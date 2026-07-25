"""Primary-source fetchers with an on-disk cache.

Every number this package reports traces back to one of three public sources:

  * UCSC hg38 REST API, ``getData/sequence``  — soft-masked reference sequence
  * UCSC hg38 REST API, ``getData/track``     — the ``gap`` and ``genomicSuperDups`` tracks
  * GIAB HG001 v4.2.1 benchmark BED           — high-confidence intervals

Responses are cached verbatim under ``cache/`` so a re-run is offline and free.
The cache is deliberately not committed: it is large, re-fetchable, and the
committed ``results/`` TSVs are what the audit trail actually rests on.

Both sources are public domain, satisfying the repository's "public data only" rule.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable
from pathlib import Path

UCSC_API = "https://api.genome.ucsc.edu"
GIAB_BENCHMARK_BED = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
)

GENOME = "hg38"
CACHE_DIR = Path(__file__).resolve().parent.parent / "cache"

# UCSC asks for courtesy on the public API endpoint. The chr20 window scan issues
# ~190 sequence requests, so pace them rather than hammering.
REQUEST_DELAY_S = 0.34
MAX_RETRIES = 4

# HTTP statuses worth retrying: request-timeout, rate-limit, and the transient
# 5xx/gateway failures. Any other status (404, 400, 403, 501…) is the request's
# own fault and fails identically on retry, so we raise at once rather than burn
# the full backoff schedule first. Retry-After (sent on some 429/503 responses)
# is deliberately not honoured — the fixed backoff below is enough for a scan
# this small; respecting it would be a separate change.
RETRYABLE_STATUS = frozenset({408, 429, 500, 502, 503, 504})

# Identifies the tool in UCSC/NCBI server logs. Anonymous ``Python-urllib`` traffic
# is what rate limiters null-route first. Deliberately carries no contact address:
# NCBI's guidelines suggest one, but that is personal data and this scan is far too
# small to warrant it.
USER_AGENT = "nf-varcall-cloud/fixture-selection"


class SourceError(RuntimeError):
    """A primary source could not be read, or returned something unexpected."""


def _fetch(url: str) -> bytes:
    """GET ``url`` and return the raw body.

    Transient failures — connection errors, timeouts, and the retryable HTTP
    statuses in :data:`RETRYABLE_STATUS` (5xx, 408, 429) — are retried with
    exponential backoff. A non-retryable HTTP status (404, 400, 403, …) is the
    request's own fault and would fail identically on retry, so it raises at
    once instead of exhausting the backoff schedule first.
    """
    last: Exception | None = None
    for attempt in range(MAX_RETRIES):
        if attempt:
            time.sleep(2**attempt)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as resp:
                body = resp.read()
            time.sleep(REQUEST_DELAY_S)
            return body
        except urllib.error.HTTPError as exc:  # subclass of URLError — catch first
            if exc.code not in RETRYABLE_STATUS:
                raise SourceError(f"HTTP {exc.code} for {url}: {exc.reason}") from exc
            last = exc
        except (urllib.error.URLError, OSError) as exc:
            last = exc
    raise SourceError(f"failed after {MAX_RETRIES} attempts: {url}\n  {last}")


def _atomic_write(path: Path, body: bytes) -> None:
    """Write ``body`` to ``path`` via a temp file + atomic rename.

    A bare ``write_bytes`` interrupted mid-write leaves a truncated file that
    still ``exists()``; a later offline run would read it as though complete.
    The rename is atomic on POSIX, so a reader sees either the old file or the
    whole new one, never a half-written one.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(body)
    tmp.replace(path)


def _cached(
    key: str, url: str, validate: Callable[[bytes], None] | None = None
) -> bytes:
    """Return the body for ``url``, reading through ``cache/<key>``.

    ``validate``, if given, must raise :class:`SourceError` for a body that must
    not be trusted. It runs *before* anything is committed to the cache, so a
    malformed HTTP 200 — a UCSC ``{"error": ...}`` payload, an FTP error page —
    is never persisted and cannot poison later offline runs. It runs again on a
    cache hit, so a file left truncated by an interrupted earlier run is dropped
    and re-fetched rather than parsed as though it were complete.
    """
    path = CACHE_DIR / key
    if path.exists():
        body = path.read_bytes()
        if validate is None:
            return body
        try:
            validate(body)
            return body
        except SourceError:
            path.unlink()  # poisoned or truncated — fall through and re-fetch
    body = _fetch(url)
    if validate is not None:
        validate(body)  # a bad fresh fetch raises here; nothing is written
    _atomic_write(path, body)
    return body


def _validate_ucsc(raw: bytes, url: str) -> None:
    """A UCSC response worth caching: parseable JSON that carries no ``error``."""
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise SourceError(f"UCSC returned non-JSON for {url}: {exc}") from None
    if "error" in payload:
        raise SourceError(f"UCSC error for {url}: {payload['error']}")


def _ucsc(endpoint: str, cache_key: str, **params: object) -> dict:
    """Call a UCSC REST endpoint. ``cache_key`` is explicit, never derived from
    parameter ordering."""
    url = f"{UCSC_API}/{endpoint}?" + urllib.parse.urlencode(
        {"genome": GENOME, **params}
    )
    raw = _cached(
        f"ucsc/{cache_key}.json", url, validate=lambda body: _validate_ucsc(body, url)
    )
    return json.loads(raw)  # already validated as parseable and error-free


def sequence(chrom: str, start: int, end: int) -> str:
    """Soft-masked reference sequence for ``chrom:start-end`` (0-based, half-open).

    Case is preserved and load-bearing: lowercase marks RepeatMasker soft-masking,
    which is how the repeat fraction is derived downstream.
    """
    payload = _ucsc(
        "getData/sequence",
        f"seq_{GENOME}_{chrom}_{start}_{end}",
        chrom=chrom,
        start=start,
        end=end,
    )
    if "dna" not in payload:
        raise SourceError(
            f"no 'dna' field in sequence response for {chrom}:{start}-{end}; "
            f"got keys {sorted(payload)}"
        )
    dna = payload["dna"]
    if len(dna) != end - start:
        raise SourceError(
            f"{chrom}:{start}-{end}: expected {end - start} bp, got {len(dna)}"
        )
    return dna


def track(name: str, chrom: str) -> list[dict]:
    """Raw feature dicts from a UCSC track, restricted to one chromosome.

    UCSC keys the payload by track name, and returns either a flat list of
    features or a per-chromosome mapping depending on the track. Both shapes are
    handled explicitly; anything else raises rather than silently yielding an
    empty list — an empty track would quietly disable an exclusion filter and
    inflate the surviving-window count.
    """
    payload = _ucsc(
        "getData/track", f"track_{GENOME}_{name}_{chrom}", track=name, chrom=chrom
    )
    if name not in payload:
        raise SourceError(
            f"track '{name}' missing from response for {chrom}; "
            f"got keys {sorted(payload)}"
        )
    data = payload[name]
    if isinstance(data, dict):  # per-chromosome mapping
        data = data.get(chrom, [])
    if not isinstance(data, list):
        raise SourceError(
            f"track '{name}' on {chrom}: expected a list of features, "
            f"got {type(data).__name__}"
        )
    if not data:
        raise SourceError(
            f"track '{name}' on {chrom} came back empty — refusing to treat this "
            f"as 'no features to exclude'"
        )
    return data


def track_intervals(name: str, chrom: str) -> list[tuple[int, int]]:
    """``(start, end)`` pairs from a UCSC track, sorted and 0-based half-open."""
    out = []
    for feature in track(name, chrom):
        if feature.get("chrom", chrom) != chrom:
            continue
        try:
            out.append((int(feature["chromStart"]), int(feature["chromEnd"])))
        except (KeyError, TypeError, ValueError) as exc:
            raise SourceError(
                f"track '{name}' on {chrom}: unparseable feature {feature!r}"
            ) from exc
    return sorted(out)


def _validate_bed(raw: bytes) -> None:
    """A body worth caching as BED: non-empty, with a tab-delimited first data
    line whose start/end are integers. Rejects an HTML error page and an empty
    or truncated download before either reaches the cache."""
    for line in raw.splitlines():
        if not line or line.startswith((b"#", b"track", b"browser")):
            continue
        cols = line.split(b"\t")
        if len(cols) < 3 or not (cols[1].isdigit() and cols[2].isdigit()):
            raise SourceError(
                f"GIAB download is not BED (first data line: {line[:80]!r})"
            )
        return
    raise SourceError("GIAB download is empty")


def _giab_bed() -> Path:
    """The genome-wide GIAB benchmark BED, downloaded once (validated, atomic)."""
    key = "giab/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
    _cached(key, GIAB_BENCHMARK_BED, validate=_validate_bed)
    return CACHE_DIR / key


def giab_highconf(chrom: str) -> list[tuple[int, int]]:
    """GIAB HG001 v4.2.1 high-confidence intervals on ``chrom``.

    Derived from the single cached genome-wide download, so asking for several
    chromosomes costs one fetch, not one per chromosome.
    """
    subset = CACHE_DIR / "giab" / f"HG001_v4.2.1_{chrom}.bed"
    if not subset.exists():
        wanted = chrom.encode()
        with _giab_bed().open("rb") as fh:
            lines = [ln for ln in fh if ln.split(b"\t", 1)[0] == wanted]
        if not lines:
            raise SourceError(f"no {chrom} intervals in the GIAB benchmark BED")
        _atomic_write(subset, b"".join(lines))
    return sorted(
        (int(p[1]), int(p[2]))
        for p in (ln.split("\t") for ln in subset.read_text().splitlines() if ln)
    )


def highconf_overlapping(chrom: str, start: int, end: int) -> list[tuple[int, int]]:
    """High-confidence intervals *overlapping* ``chrom:start-end``, untruncated.

    Intervals are returned at their full extent, so the first and last may reach
    outside the window. That matters: the site filters test membership against
    the real interval bounds, and truncating here would fabricate high-confidence
    edges exactly where the window boundary falls.
    """
    return [(s, e) for s, e in giab_highconf(chrom) if s < end and e > start]


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()

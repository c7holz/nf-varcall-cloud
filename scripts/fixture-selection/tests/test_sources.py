"""Offline unit tests for :mod:`fixtures.sources`.

Run from ``scripts/fixture-selection/``::

    python3 -m unittest discover -v

Stdlib only, no venv, and never in CI — same rule as the package under test.
Every test here is offline: the network boundary is stubbed, and
:class:`CacheTestCase` additionally makes a real request *fail loudly* rather
than succeed quietly. What the live sources actually return is not this file's
concern; that is the ``reproduce`` command's job, which re-derives the committed
``results/`` end to end.

These tests cover the places where a silent bug corrupts the audit — validation
that must reject a bad body before it reaches the cache, interval parsing, and
the retry classification in ``_fetch``.

Conventions
-----------

* **URLs are not regex.** Source URLs contain ``?`` and ``&``, which
  ``assertRaisesRegex`` would interpret. Assertions on a message carrying a URL
  or path use ``assertRaises(...) as cm`` plus ``assertIn``; ``assertRaisesRegex``
  is reserved for regex-safe literal fragments.
* **``subTest`` loops that assert call counts reset their mocks** at the top of
  each iteration, otherwise counts accumulate across cases and the assertion
  means nothing.
* **urllib capitalizes header lookups.** ``Request.get_header("User-Agent")``
  returns ``None``; only ``get_header("User-agent")`` finds it — a test using
  the obvious spelling would pass against nothing.
* **Pass an explicit no-op validator** (``lambda body: None``) in the ``_cached``
  tests that are not exercising validation.
"""

from __future__ import annotations

import io
import unittest
import urllib.error
import urllib.parse
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

from fixtures import sources

# A stand-in source URL. Carries a query string on purpose: '?' and '&' are
# regex metacharacters, which is why assertions on messages containing a URL use
# assertRaises + assertIn rather than assertRaisesRegex.
URL = "https://api.genome.ucsc.edu/getData/sequence?genome=hg38&chrom=chr20"

# Reason phrase for every stubbed HTTPError. Deliberately not a real phrase
# ("Not Found", "Service Unavailable"): if this string ever surfaces from code
# under test, it is unambiguously the fake talking and not a live response.
FAKE_REASON = "stubbed HTTP error (unit test)"


def http_error(code: int) -> urllib.error.HTTPError:
    """An ``HTTPError`` carrying ``code``, with ``.reason`` reading back
    :data:`FAKE_REASON`."""
    # Both odd details below are load-bearing. HTTPError doubles as a response
    # object, so an unclosed one emits a ResourceWarning mid-run, blamed on
    # whichever test is executing; closing it keeps code and reason, which is
    # all _fetch reads. The explicit buffer is what makes close() safe before
    # Python 3.13, which is the first version to substitute one for fp=None.
    exc = urllib.error.HTTPError(URL, code, FAKE_REASON, {}, io.BytesIO())
    exc.close()
    return exc


class FakeResponse:
    """Minimal stand-in for what ``urlopen`` returns: a context manager with
    ``read()``. That is the whole surface ``sources._fetch`` uses."""

    def __init__(self, body: bytes) -> None:
        self._body = body

    def __enter__(self) -> FakeResponse:
        return self

    def __exit__(self, *exc: object) -> bool:
        return False

    def read(self) -> bytes:
        return self._body


def query_of(url: str) -> dict[str, list[str]]:
    """Parsed query parameters of ``url``, so URL assertions do not depend on
    the ordering ``urlencode`` happens to produce."""
    return urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)


class CacheTestCase(unittest.TestCase):
    """Base for tests that touch the cache, with the network taken away.

    ``setUp`` points ``sources.CACHE_DIR`` at a fresh temporary directory —
    which works because ``_cached`` and ``giab_highconf`` read that global at
    call time — and replaces ``urllib.request.urlopen`` with a land mine, so a
    code path that slips past a stubbed ``_fetch`` fails instead of quietly
    reaching UCSC. Subclasses that stub ``_fetch`` call ``super().setUp()``
    first, then add their own patch.
    """

    def setUp(self) -> None:
        tmp = TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.cache_dir = Path(tmp.name)

        cache_dir = mock.patch.object(sources, "CACHE_DIR", self.cache_dir)
        cache_dir.start()
        self.addCleanup(cache_dir.stop)

        urlopen = mock.patch("urllib.request.urlopen")
        self.urlopen = urlopen.start()
        self.addCleanup(urlopen.stop)
        self.urlopen.side_effect = AssertionError("network access in a unit test")


class ValidateUcscTests(unittest.TestCase):
    """``_validate_ucsc`` — the gate a UCSC body passes before it is cached."""

    def test_accepts_valid_payload(self) -> None:
        self.assertIsNone(sources._validate_ucsc(b'{"dna": "ACGT"}', URL))

    def test_rejects_non_json(self) -> None:
        with self.assertRaises(sources.SourceError) as cm:
            sources._validate_ucsc(b"<html>503 Service Unavailable</html>", URL)
        self.assertIn(URL, str(cm.exception))

    def test_rejects_error_payload_with_http_200(self) -> None:
        # The case the validator exists for: UCSC reports failure in a 200 body,
        # so _fetch sees a perfectly good response and would cache the error.
        with self.assertRaisesRegex(sources.SourceError, "no such table"):
            sources._validate_ucsc(b'{"error": "no such table"}', URL)


class ValidateBedTests(unittest.TestCase):
    """``_validate_bed`` — the gate a GIAB body passes before it is cached."""

    def test_accepts_minimal_bed(self) -> None:
        self.assertIsNone(sources._validate_bed(b"chr20\t100\t200\n"))

    def test_accepts_extra_columns(self) -> None:
        self.assertIsNone(sources._validate_bed(b"chr20\t100\t200\tname\t0\t+\n"))

    def test_skips_header_lines(self) -> None:
        raw = (
            b"# a comment\n"
            b'track name="HG001"\n'
            b"browser position chr20\n"
            b"chr20\t100\t200\n"
        )
        self.assertIsNone(sources._validate_bed(raw))

    def test_rejects_html_error_page(self) -> None:
        # The case the validator exists for: NCBI's FTP-over-HTTP serves error
        # pages with status 200, so _fetch sees a perfectly good response.
        with self.assertRaisesRegex(sources.SourceError, "not BED"):
            sources._validate_bed(b"<!DOCTYPE html>\n<html>404</html>\n")

    def test_rejects_empty_download(self) -> None:
        for raw in (b"", b"\n\n", b"# only a comment\n"):
            with self.subTest(raw=raw):
                with self.assertRaisesRegex(sources.SourceError, "empty"):
                    sources._validate_bed(raw)

    def test_rejects_non_integer_coordinates(self) -> None:
        for raw in (b"chr20\t100\n", b"chr20\tfoo\t200\n", b"chr20\t100\tbar\n"):
            with self.subTest(raw=raw):
                with self.assertRaisesRegex(sources.SourceError, "not BED"):
                    sources._validate_bed(raw)

    def test_only_the_first_data_line_is_checked(self) -> None:
        # Deliberate: the guard's reach ends at the first data line, so a body
        # truncated mid-file passes here. Detecting that is a digest's job, not
        # this validator's — see the phase-4 per-source sha256 requirement.
        self.assertIsNone(sources._validate_bed(b"chr20\t100\t200\nnot a bed line\n"))

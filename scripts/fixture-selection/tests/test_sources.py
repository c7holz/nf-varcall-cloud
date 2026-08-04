"""Offline unit tests for :mod:`fixtures.sources`.

Run from ``scripts/fixture-selection/``::

    python3 -m unittest discover -v

Stdlib only, no venv, never in CI — the same rule as the package under test.
Everything here is offline: the network boundary is stubbed, and
:class:`CacheTestCase` makes a stray real request *fail loudly* rather than
succeed quietly. What the live sources actually return is the ``reproduce``
command's concern, not this file's. A test earns its place here by covering a
spot where a silent bug would corrupt the audit.

Conventions
-----------

* **URLs are not regex.** They carry ``?`` and ``&``, which
  ``assertRaisesRegex`` would interpret. Assert on a message containing a URL
  or a path with ``assertRaises(...) as cm`` plus ``assertIn``.
* **Reset the mocks inside a ``subTest`` loop that counts calls**, or the
  counts accumulate across cases and the assertion means nothing.
* **Pass an explicit no-op validator** (``lambda body: None``) to ``_cached``
  in tests that are not exercising validation.
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

    def test_rejects_non_object_json(self) -> None:
        # _cached drops a cache entry only when the validator raises SourceError,
        # so a TypeError here wedges the entry instead of healing it.
        for raw in (b"0", b"null", b'"error"', b"[1, 2]"):
            with self.subTest(raw=raw):
                with self.assertRaises(sources.SourceError) as cm:
                    sources._validate_ucsc(raw, URL)
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


class FetchTests(unittest.TestCase):
    """``_fetch`` — the single network boundary, and how it classifies failure.

    Deliberately not a :class:`CacheTestCase` subclass: this is the one place
    that *should* call ``urlopen``, so it owns a real mock rather than
    inheriting the land mine. ``time.sleep`` is patched because ``_fetch``
    otherwise waits 2+4+8 s across a full retry schedule.
    """

    def setUp(self) -> None:
        urlopen = mock.patch("urllib.request.urlopen")
        self.urlopen = urlopen.start()
        self.addCleanup(urlopen.stop)

        sleep = mock.patch("time.sleep")
        self.sleep = sleep.start()
        self.addCleanup(sleep.stop)

    def test_returns_body_on_first_try(self) -> None:
        self.urlopen.return_value = FakeResponse(b"payload")
        self.assertEqual(sources._fetch(URL), b"payload")
        self.assertEqual(self.urlopen.call_count, 1)

    def test_sends_project_user_agent(self) -> None:
        self.urlopen.return_value = FakeResponse(b"payload")
        sources._fetch(URL)

        req = self.urlopen.call_args.args[0]
        # urllib capitalizes header names on the way in, so "User-Agent" reads
        # back as None and an assertion using it would pass against nothing.
        self.assertEqual(req.get_header("User-agent"), sources.USER_AGENT)
        # Identify the tool, carry no contact address: anonymous urllib gets
        # rate-limited first, but an email in the header is personal data.
        self.assertNotIn("@", sources.USER_AGENT)

    def test_fast_fails_on_non_retryable_status(self) -> None:
        # The client's own fault — retrying cannot change the answer. Also
        # guards the except-arm ordering: HTTPError ⊂ URLError ⊂ OSError, so
        # the HTTPError arm must stay first or every 404 starts retrying.
        for code in (400, 403, 404, 501):
            with self.subTest(code=code):
                self.urlopen.reset_mock()
                self.sleep.reset_mock()
                self.urlopen.side_effect = http_error(code)

                with self.assertRaisesRegex(sources.SourceError, f"HTTP {code}"):
                    sources._fetch(URL)

                self.assertEqual(self.urlopen.call_count, 1)
                self.assertEqual(self.sleep.call_count, 0)

    def test_retries_retryable_status_then_succeeds(self) -> None:
        # Includes 408 and 429 — the two 4xx codes that must not take the
        # fast-fail branch above.
        for code in sorted(sources.RETRYABLE_STATUS):
            with self.subTest(code=code):
                self.urlopen.reset_mock()
                self.urlopen.side_effect = [
                    http_error(code),
                    http_error(code),
                    FakeResponse(b"ok"),
                ]

                self.assertEqual(sources._fetch(URL), b"ok")
                self.assertEqual(self.urlopen.call_count, 3)

    def test_retries_transient_connection_errors(self) -> None:
        self.urlopen.side_effect = [
            urllib.error.URLError("timed out"),
            FakeResponse(b"ok"),
        ]

        self.assertEqual(sources._fetch(URL), b"ok")
        self.assertEqual(self.urlopen.call_count, 2)

    def test_gives_up_after_max_retries(self) -> None:
        self.urlopen.side_effect = urllib.error.URLError("timed out")

        with self.assertRaises(sources.SourceError) as cm:
            sources._fetch(URL)

        message = str(cm.exception)
        self.assertIn(f"failed after {sources.MAX_RETRIES} attempts", message)
        self.assertIn(URL, message)
        # The last error's own text survives into the message: without it the
        # caller learns that every attempt failed but never why.
        self.assertIn("timed out", message)
        self.assertEqual(self.urlopen.call_count, sources.MAX_RETRIES)

    def test_backoff_is_exponential_and_paced(self) -> None:
        self.urlopen.side_effect = [
            urllib.error.URLError("timed out"),
            urllib.error.URLError("timed out"),
            FakeResponse(b"ok"),
        ]

        sources._fetch(URL)

        # Two backoffs before the retries, then the courtesy delay UCSC is owed
        # once the request finally succeeds.
        self.assertEqual(
            self.sleep.mock_calls,
            [mock.call(2), mock.call(4), mock.call(sources.REQUEST_DELAY_S)],
        )


class CachedTests(CacheTestCase):
    """``_cached`` — the read-through cache, and what it refuses to keep.

    Validation runs on both paths on purpose, and these tests pin both: a body
    that fails the gate must never end up on disk, because the cache is what a
    later offline run trusts in place of the source.
    """

    KEY = "ucsc/x.json"

    def setUp(self) -> None:
        super().setUp()
        fetch = mock.patch.object(sources, "_fetch")
        self.fetch = fetch.start()
        self.addCleanup(fetch.stop)

    def prewrite(self, body: bytes) -> Path:
        """Seed the cache entry as an earlier run would have left it."""
        path = self.cache_dir / self.KEY
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
        return path

    def test_fetch_on_miss_writes_cache(self) -> None:
        self.fetch.return_value = b"body"

        body = sources._cached(self.KEY, URL, validate=lambda body: None)

        self.assertEqual(body, b"body")
        self.assertEqual((self.cache_dir / self.KEY).read_bytes(), b"body")
        self.fetch.assert_called_once_with(URL)

    def test_creates_nested_parent_directories(self) -> None:
        # The key carries a subdirectory that no earlier run has created.
        nested = self.cache_dir / "ucsc"
        self.assertFalse(nested.exists())
        self.fetch.return_value = b"body"

        sources._cached(self.KEY, URL, validate=lambda body: None)

        self.assertTrue(nested.is_dir())
        # _atomic_write renames its temp file rather than leaving it behind; a
        # stray .tmp would be read by nothing but would mean the rename is gone.
        self.assertEqual(list(nested.glob("*.tmp")), [])

    def test_cache_hit_does_not_fetch(self) -> None:
        self.prewrite(b"body")

        body = sources._cached(self.KEY, URL, validate=lambda body: None)

        self.assertEqual(body, b"body")
        self.fetch.assert_not_called()

    def test_validates_on_cache_hit(self) -> None:
        self.prewrite(b"body")
        validate = mock.Mock(return_value=None)

        body = sources._cached(self.KEY, URL, validate=validate)

        self.assertEqual(body, b"body")
        validate.assert_called_once_with(b"body")
        self.fetch.assert_not_called()

    def test_drops_poisoned_cache_entry_and_refetches(self) -> None:
        # An error page an earlier run cached before the validator existed, or
        # a file left truncated by an interrupted run: rejecting it on read is
        # what stops one bad download from outliving the run that made it.
        path = self.prewrite(b"<html>")
        self.fetch.return_value = b"good"

        def validate(body: bytes) -> None:
            if body == b"<html>":
                raise sources.SourceError("poisoned")

        body = sources._cached(self.KEY, URL, validate=validate)

        self.assertEqual(body, b"good")
        self.assertEqual(path.read_bytes(), b"good")
        self.fetch.assert_called_once_with(URL)

    def test_does_not_cache_a_bad_fresh_fetch(self) -> None:
        self.fetch.return_value = b"<html>"

        def validate(body: bytes) -> None:
            raise sources.SourceError("not JSON")

        with self.assertRaisesRegex(sources.SourceError, "not JSON"):
            sources._cached(self.KEY, URL, validate=validate)

        # Nothing written at all: a raise the caller retries past would
        # otherwise leave the bad body to be served offline as though good.
        self.assertFalse((self.cache_dir / self.KEY).exists())
        self.assertEqual(list(self.cache_dir.rglob("*.tmp")), [])

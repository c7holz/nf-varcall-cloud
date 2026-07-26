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
    """An ``HTTPError`` carrying ``code``. ``fp=None`` is supported, and
    ``.reason`` reads back :data:`FAKE_REASON`."""
    return urllib.error.HTTPError(URL, code, FAKE_REASON, {}, None)


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

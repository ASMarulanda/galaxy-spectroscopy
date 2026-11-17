"""Retrieval: identifiers, URLs, and every way a download can go wrong.

No test here touches the network. The archive is replaced by ``FakeOpener``, which
answers the same two calls ``UrllibOpener`` makes (``head`` and ``get``) and can be
told to fail in the specific ways the real one fails: a 404 for an object that does
not exist, a 503 for a server having a bad afternoon, a socket timeout, a body that
stops early, and a server that ignores ``Range``.

``sleep`` is injected everywhere so backoff is exercised without waiting.
"""

from __future__ import annotations

import hashlib
import io
import urllib.error
from email.message import Message

import pytest

import download_sdss as D

BODY = bytes(range(256)) * 4            # 1024 deterministic bytes
SHA = hashlib.sha256(BODY).hexdigest()


class Resp:
    def __init__(self, body: bytes, status: int = 200):
        self._buf = io.BytesIO(body)
        self.status = status
        self.closed = False

    def read(self, n: int = -1) -> bytes:
        return self._buf.read(n)

    def close(self):
        self.closed = True


class FakeOpener:
    """A scripted stand-in for UrllibOpener."""

    def __init__(self, body: bytes = BODY, *, head_status: int = 200,
                 content_length: int | None = -1, accept_ranges: str | None = "bytes",
                 head_failures: int = 0, get_failures: int = 0,
                 ignore_range: bool = False, truncate_to: int | None = None):
        self.body = body
        self.head_status = head_status
        self.content_length = len(body) if content_length == -1 else content_length
        self.accept_ranges = accept_ranges
        self.head_failures = head_failures
        self.get_failures = get_failures
        self.ignore_range = ignore_range
        self.truncate_to = truncate_to
        self.head_calls: list[str] = []
        self.get_calls: list[int] = []

    def head(self, url: str) -> D.HeadInfo:
        self.head_calls.append(url)
        if self.head_failures:
            self.head_failures -= 1
            raise urllib.error.HTTPError(url, 503, "unavailable", Message(), None)
        return D.HeadInfo(status=self.head_status, content_length=self.content_length,
                          etag='"abc"', last_modified="Tue, 01 Jan 2013 00:00:00 GMT",
                          accept_ranges=self.accept_ranges)

    def get(self, url: str, start: int = 0):
        self.get_calls.append(start)
        if self.get_failures:
            self.get_failures -= 1
            raise urllib.error.URLError("connection reset")
        if self.ignore_range:
            return Resp(self.body[: self.truncate_to] if self.truncate_to else self.body, 200)
        if self.truncate_to is not None:
            # a connection that dies: the first response carries the first N bytes,
            # every later read returns EOF, exactly as a dead socket does
            if start >= self.truncate_to:
                return Resp(b"", 206 if start else 200)
            return Resp(self.body[start: self.truncate_to], 206 if start else 200)
        return Resp(self.body[start:], 206 if start else 200)


NO_SLEEP = {"sleep": lambda _s: None}


# ------------------------------------------------------------------ identifiers
@pytest.mark.parametrize("value", [36, "36", " 36 ", "+36"])
def test_integer_identifiers_are_accepted(value):
    assert D.validate_identifier(value, "fiber", (1, 1000)) == 36


@pytest.mark.parametrize("value", ["abc", "", " ", "36; rm -rf /", "3 6", "36.0",
                                   "0x24", "-5", "1e3", None, True, False])
def test_non_integer_identifiers_are_rejected(value):
    with pytest.raises(D.IdentifierError):
        D.validate_identifier(value, "fiber", (1, 1000))


@pytest.mark.parametrize("plate,mjd,fiber", [(0, 51986, 36), (10000, 51986, 36),
                                             (294, 1, 36), (294, 70000, 36),
                                             (294, 51986, 0), (294, 51986, 1001)])
def test_out_of_range_identifiers_are_rejected(plate, mjd, fiber):
    with pytest.raises(D.IdentifierError, match="outside the documented range"):
        D.validate_triple(plate, mjd, fiber)


def test_url_and_filename_construction():
    assert D.url_for(294, 51986, 36) == (
        "https://data.sdss.org/sas/dr18/prior-surveys/sdss2-dr8-sdss/spectro/"
        "redux/26/spectra/0294/spec-0294-51986-0036.fits")
    assert D.filename_for(294, 51986, 36) == "spec-0294-51986-0036.fits"
    # both the directory and the file name are zero-padded
    assert D.url_for(7, 51234, 5).endswith("/0007/spec-0007-51234-0005.fits")
    assert D.case_id_for(294, 51986, 36) == "SDSS-0294-51986-0036"


def test_url_is_built_only_from_validated_integers():
    """A hostile string cannot reach the URL: it fails validation first."""
    with pytest.raises(D.IdentifierError):
        D.url_for(294, 51986, "36&x=../../etc/passwd")     # type: ignore[arg-type]


# ------------------------------------------------------------------------ http
def test_missing_object_is_404_and_is_not_retried(tmp_path):
    opener = FakeOpener(head_status=404)
    with pytest.raises(D.FetchError) as exc:
        D.fetch_one(294, 51986, 36, tmp_path, opener=opener, **NO_SLEEP)
    assert exc.value.kind == "notfound"
    assert exc.value.exit_code == 1
    assert len(opener.head_calls) == 1          # no pointless repetition
    assert not list(tmp_path.glob("*.fits"))


def test_server_errors_are_retried_then_succeed(tmp_path):
    opener = FakeOpener(head_failures=2)
    result = D.fetch_one(294, 51986, 36, tmp_path, opener=opener, **NO_SLEEP)
    assert len(opener.head_calls) == 3
    assert result.status == "fetched"
    assert result.sha256 == SHA


def test_retries_are_bounded_and_report_network(tmp_path):
    opener = FakeOpener(head_failures=99)
    with pytest.raises(D.FetchError) as exc:
        D.fetch_one(294, 51986, 36, tmp_path, opener=opener, attempts=3, **NO_SLEEP)
    assert exc.value.kind == "network"
    assert exc.value.exit_code == 3
    assert len(opener.head_calls) == 3


def test_transport_failure_during_body_is_retried(tmp_path):
    opener = FakeOpener(get_failures=1)
    result = D.fetch_one(294, 51986, 36, tmp_path, opener=opener, **NO_SLEEP)
    assert result.status == "fetched"
    assert (tmp_path / "spec-0294-51986-0036.fits").read_bytes() == BODY


def test_truncated_download_is_detected_and_the_part_is_kept(tmp_path):
    opener = FakeOpener(truncate_to=700)        # content-length says 1024
    with pytest.raises(D.FetchError) as exc:
        D.fetch_one(294, 51986, 36, tmp_path, opener=opener, **NO_SLEEP)
    assert exc.value.kind == "verify"
    assert "incomplete" in str(exc.value)
    assert (tmp_path / "spec-0294-51986-0036.fits.part").exists()
    assert not (tmp_path / "spec-0294-51986-0036.fits").exists()


def test_sha256_mismatch_fails_and_deletes_the_part(tmp_path):
    opener = FakeOpener()
    with pytest.raises(D.FetchError) as exc:
        D.fetch_one(294, 51986, 36, tmp_path, opener=opener,
                    expected_sha256="0" * 64, **NO_SLEEP)
    assert exc.value.kind == "verify"
    assert not (tmp_path / "spec-0294-51986-0036.fits.part").exists()
    assert not (tmp_path / "spec-0294-51986-0036.fits").exists()


def test_sha256_that_matches_is_accepted(tmp_path):
    result = D.fetch_one(294, 51986, 36, tmp_path, opener=FakeOpener(),
                         expected_sha256=SHA, **NO_SLEEP)
    assert result.sha256 == SHA


def test_resume_continues_from_the_part_file(tmp_path):
    part = tmp_path / "spec-0294-51986-0036.fits.part"
    part.write_bytes(BODY[:400])
    opener = FakeOpener()
    result = D.fetch_one(294, 51986, 36, tmp_path, opener=opener, **NO_SLEEP)
    assert opener.get_calls == [400]            # asked for the remainder only
    assert result.status == "resumed"
    assert (tmp_path / "spec-0294-51986-0036.fits").read_bytes() == BODY
    assert result.http["resumed_from"] == 400


def test_resume_restarts_when_the_server_ignores_range(tmp_path):
    """A server that advertises Range and then ignores it must not strand a .part."""
    part = tmp_path / "spec-0294-51986-0036.fits.part"
    part.write_bytes(BODY[:400])
    opener = FakeOpener(ignore_range=True)
    result = D.fetch_one(294, 51986, 36, tmp_path, opener=opener, **NO_SLEEP)
    assert result.status == "fetched"
    assert result.http["resumed_from"] == 0          # the second attempt started over
    assert (tmp_path / "spec-0294-51986-0036.fits").read_bytes() == BODY


def test_no_resume_restarts_from_zero(tmp_path):
    (tmp_path / "spec-0294-51986-0036.fits.part").write_bytes(BODY[:400])
    opener = FakeOpener()
    result = D.fetch_one(294, 51986, 36, tmp_path, opener=opener, resume=False, **NO_SLEEP)
    assert opener.get_calls == [0]
    assert (tmp_path / "spec-0294-51986-0036.fits").read_bytes() == BODY


def test_part_file_is_gone_after_success(tmp_path):
    D.fetch_one(294, 51986, 36, tmp_path, opener=FakeOpener(), **NO_SLEEP)
    assert not list(tmp_path.glob("*.part"))


# ---------------------------------------------------------------- no overwriting
def test_existing_file_is_not_replaced(tmp_path):
    dest = tmp_path / "spec-0294-51986-0036.fits"
    dest.write_bytes(b"an older copy")
    opener = FakeOpener()
    result = D.fetch_one(294, 51986, 36, tmp_path, opener=opener, **NO_SLEEP)
    assert result.status == "present"
    assert opener.head_calls == [] and opener.get_calls == []
    assert dest.read_bytes() == b"an older copy"


def test_existing_file_with_wrong_hash_is_refused(tmp_path):
    dest = tmp_path / "spec-0294-51986-0036.fits"
    dest.write_bytes(b"corrupt")
    with pytest.raises(D.FetchError, match="refusing to overwrite"):
        D.fetch_one(294, 51986, 36, tmp_path, opener=FakeOpener(),
                    expected_sha256=SHA, **NO_SLEEP)
    assert dest.read_bytes() == b"corrupt"


def test_force_replaces_an_existing_file(tmp_path):
    dest = tmp_path / "spec-0294-51986-0036.fits"
    dest.write_bytes(b"old")
    result = D.fetch_one(294, 51986, 36, tmp_path, opener=FakeOpener(), force=True, **NO_SLEEP)
    assert result.status == "fetched"
    assert dest.read_bytes() == BODY


# ------------------------------------------------------------------- list input
def test_list_csv_accepts_the_catalogue_column_names(tmp_path):
    path = tmp_path / "objects.csv"
    path.write_text("plate,fiberID,mjd\n294,36,51986\n408,446,51821\n")
    assert D.read_list(path) == [(294, 51986, 36), (408, 51821, 446)]


def test_list_csv_reports_the_offending_line(tmp_path):
    path = tmp_path / "objects.csv"
    path.write_text("plate,mjd,fiber\n294,51986,36\n294,51986,oops\n")
    with pytest.raises(D.FetchError) as exc:
        D.read_list(path)
    assert exc.value.kind == "usage"
    assert "line 3" in str(exc.value)


def test_list_csv_needs_the_right_columns(tmp_path):
    path = tmp_path / "objects.csv"
    path.write_text("a,b,c\n1,2,3\n")
    with pytest.raises(D.FetchError, match="needs a column"):
        D.read_list(path)


# ---------------------------------------------------------------- exit codes
@pytest.mark.parametrize("kind,code", [("usage", 1), ("notfound", 1), ("verify", 2),
                                       ("network", 3), ("fits", 4)])
def test_exit_codes_match_the_documented_table(kind, code):
    assert D.FetchError("x", kind).exit_code == code


def test_cli_rejects_a_bad_identifier_before_any_request(capsys):
    code = D.main(["--plate", "294", "--mjd", "51986", "--fiber", "abc",
                   "--out", "/tmp/should-not-be-used"])
    assert code == 1
    assert "plain integer" in capsys.readouterr().err


def test_cli_fetch_then_convert_without_network(tmp_path, monkeypatch, capsys):
    """The whole CLI path, with the archive replaced by the fake."""
    monkeypatch.setattr(D, "UrllibOpener", lambda **kw: FakeOpener(body=BODY))
    code = D.main(["--plate", "294", "--mjd", "51986", "--fiber", "36",
                   "--out", str(tmp_path)])
    # BODY is not a FITS file, so the contract check must refuse it: exit 4,
    # and the failure must be recorded rather than silently skipped.
    assert code == 4
    manifest = (tmp_path / "provenance" / "fetch_manifest.json").read_text()
    assert '"failed": 1' in manifest
    assert "not a readable FITS" in capsys.readouterr().err

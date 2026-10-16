"""Focused regression tests for the build-time commit hash helper.

``get_commit_hash`` runs from ``meson.build`` on every build, and meson invokes
it with ``check: true``. Anything it raises therefore aborts the build outright,
and whatever it writes is imported by ``cosmic/__init__.py``, so both the
failure modes and the generated file are worth pinning down.
"""

import ast
import subprocess

import pytest

from cosmic import get_commit_hash as gch


VALID_SHA1 = "8de569b3f71a7f6b0450ba014be072175c5d670d"
VALID_SHA256 = "3f" * 32


def _refuse_to_run(*args, **kwargs):
    raise AssertionError("git must not be consulted when the hash is supplied")


def _fake_git(returncode, stdout=b""):
    def run(*args, **kwargs):
        return subprocess.CompletedProcess(args, returncode, stdout=stdout)

    return run


def _read_back(tmp_path):
    """Parse the generated module, asserting it holds exactly one assignment."""
    written = (tmp_path / "src" / "cosmic" / "_commit_hash.py").read_text()
    tree = ast.parse(written)
    assert len(tree.body) == 1, f"expected a single statement, got {written!r}"
    (statement,) = tree.body
    assert isinstance(statement, ast.Assign)
    assert [target.id for target in statement.targets] == ["COMMIT_HASH"]
    return ast.literal_eval(statement.value)


@pytest.fixture(autouse=True)
def _prepare_tree(tmp_path, monkeypatch):
    package_dir = tmp_path / "src" / "cosmic"
    package_dir.mkdir(parents=True)
    (tmp_path / ".git").mkdir()
    monkeypatch.setattr(gch, "__file__", str(package_dir / "get_commit_hash.py"))
    monkeypatch.chdir(tmp_path)


def test_environment_override_short_circuits_git(monkeypatch):
    monkeypatch.setenv("COSMIC_COMMIT_HASH", VALID_SHA1)
    monkeypatch.setattr(gch.subprocess, "run", _refuse_to_run)

    assert gch.get_commit_hash() == VALID_SHA1


def test_environment_override_is_stripped(monkeypatch):
    monkeypatch.setenv("COSMIC_COMMIT_HASH", f"  {VALID_SHA1}\n")
    monkeypatch.setattr(gch.subprocess, "run", _refuse_to_run)

    assert gch.get_commit_hash() == VALID_SHA1


@pytest.mark.parametrize("value", [VALID_SHA1.upper(), VALID_SHA256])
def test_uppercase_and_sha256_object_ids_are_accepted(monkeypatch, value):
    monkeypatch.setenv("COSMIC_COMMIT_HASH", value)
    monkeypatch.setattr(gch.subprocess, "run", _refuse_to_run)

    assert gch.get_commit_hash() == value


@pytest.mark.parametrize(
    "value",
    [
        "not-a-hash",
        VALID_SHA1[:-1],
        VALID_SHA1 + "0",
        "g" * 40,
        VALID_SHA1 + '"; import os',
    ],
    ids=["nonsense", "too-short", "too-long", "non-hex", "trailing-payload"],
)
def test_malformed_environment_override_is_rejected(monkeypatch, value):
    monkeypatch.setenv("COSMIC_COMMIT_HASH", value)

    with pytest.raises(ValueError):
        gch.get_commit_hash()


def test_blank_environment_override_falls_through_to_git(monkeypatch):
    monkeypatch.setenv("COSMIC_COMMIT_HASH", "   ")
    monkeypatch.setattr(gch.subprocess, "run", _fake_git(0, VALID_SHA1.encode() + b"\n"))

    assert gch.get_commit_hash() == VALID_SHA1


def test_missing_git_binary_yields_an_empty_hash(monkeypatch):
    """Container images build without git installed; this must not raise.

    Regression: an earlier revision let FileNotFoundError escape, which aborted
    the whole meson build rather than degrading to an empty hash.
    """
    monkeypatch.delenv("COSMIC_COMMIT_HASH", raising=False)

    def missing_binary(*args, **kwargs):
        raise FileNotFoundError(2, "No such file or directory: 'git'")

    monkeypatch.setattr(gch.subprocess, "run", missing_binary)

    assert gch.get_commit_hash() == ""


def test_outside_a_repository_yields_an_empty_hash(monkeypatch):
    monkeypatch.delenv("COSMIC_COMMIT_HASH", raising=False)
    monkeypatch.setattr(gch.subprocess, "run", _fake_git(128))

    assert gch.get_commit_hash() == ""


def test_unexpected_git_output_is_rejected(monkeypatch):
    monkeypatch.delenv("COSMIC_COMMIT_HASH", raising=False)
    monkeypatch.setattr(gch.subprocess, "run", _fake_git(0, b"ref: refs/heads/develop\n"))

    with pytest.raises(ValueError):
        gch.get_commit_hash()


def _init_repository(path, message="fixture"):
    path.mkdir(exist_ok=True)
    subprocess.run(["git", "init", "-q", str(path)], check=True)
    subprocess.run(
        ["git", "-C", str(path), "-c", "user.name=Test", "-c",
         "user.email=test@example.invalid", "-c", "commit.gpgsign=false",
         "commit", "--allow-empty", "-qm", message],
        check=True,
    )
    return subprocess.check_output(
        ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
    ).strip()


def test_git_uses_source_tree_from_another_repository(tmp_path, monkeypatch):
    monkeypatch.delenv("COSMIC_COMMIT_HASH", raising=False)
    expected = _init_repository(tmp_path)
    other = tmp_path / "unrelated"
    _init_repository(other, message="other")
    monkeypatch.chdir(other)

    assert gch.get_commit_hash() == expected


def test_source_archive_does_not_inherit_parent_repository(tmp_path, monkeypatch):
    monkeypatch.delenv("COSMIC_COMMIT_HASH", raising=False)
    _init_repository(tmp_path)
    archive = tmp_path / "archive" / "src" / "cosmic"
    archive.mkdir(parents=True)
    monkeypatch.setattr(gch, "__file__", str(archive / "get_commit_hash.py"))
    monkeypatch.chdir(archive)

    assert gch.get_commit_hash() == ""


@pytest.mark.parametrize("value", [VALID_SHA1, ""], ids=["hash", "empty"])
def test_written_module_round_trips(tmp_path, monkeypatch, value):
    elsewhere = tmp_path / "unrelated"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    gch.write_commit_hash_to_file(value)

    assert _read_back(tmp_path) == value


def test_written_module_contains_hostile_values(tmp_path, monkeypatch):
    """The generated file is imported, so a value must never become code.

    ``write_commit_hash_to_file`` is reachable with an unvalidated value, so it
    quotes with repr rather than interpolating into a string literal.
    """
    hostile = 'abc"\nINJECTED = True\n"'

    gch.write_commit_hash_to_file(hostile)

    assert _read_back(tmp_path) == hostile

    # cosmic/__init__.py imports this module, so execute it the way the
    # interpreter would and confirm the payload stayed inert data.
    namespace = {}
    exec((tmp_path / "src" / "cosmic" / "_commit_hash.py").read_text(), namespace)
    assert namespace["COMMIT_HASH"] == hostile
    assert "INJECTED" not in namespace

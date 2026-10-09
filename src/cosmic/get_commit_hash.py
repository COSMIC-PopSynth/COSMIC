import os
import re
import subprocess
from pathlib import Path


_COMMIT_HASH_RE = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z", re.IGNORECASE)


def _validate_commit_hash(commit_hash):
    if not _COMMIT_HASH_RE.fullmatch(commit_hash):
        raise ValueError(
            "commit hash must be a full 40- or 64-character hexadecimal Git object ID"
        )
    return commit_hash


def get_commit_hash():
    # Builds that happen outside a git checkout (container images, sdists) can
    # pass the hash in directly. Fall back to git, and to an empty string when
    # git is unavailable or this is not a repository.
    env_hash = os.environ.get('COSMIC_COMMIT_HASH', '').strip()
    if env_hash:
        return _validate_commit_hash(env_hash)

    source_root = Path(__file__).resolve().parents[2]
    # Do not mistake a source archive's enclosing repository for COSMIC.
    if not (source_root / '.git').exists():
        return ''

    try:
        result = subprocess.run(
            ['git', 'rev-parse', 'HEAD'],
            cwd=source_root,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError:
        return ''

    if result.returncode != 0:
        return ''

    return _validate_commit_hash(result.stdout.decode('utf-8').strip())


def write_commit_hash_to_file(commit_hash):
    with Path(__file__).resolve().with_name('_commit_hash.py').open('w') as f:
        f.write(f'COMMIT_HASH = {commit_hash!r}\n')


if __name__ == "__main__":
    commit_hash = get_commit_hash()
    write_commit_hash_to_file(commit_hash)

import subprocess

def get_commit_hash():
    # Run git command to get the latest commit hash
    result = subprocess.run(['git', 'rev-parse', 'HEAD'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    commit_hash = result.stdout.decode('utf-8').strip()
    return commit_hash

def write_commit_hash_to_file(commit_hash):
    with open('./src/cosmic/_commit_hash.py', 'w') as f:
        f.write(f'COMMIT_HASH = "{commit_hash}"\n')

if __name__ == "__main__":
    commit_hash = get_commit_hash()
    write_commit_hash_to_file(commit_hash)
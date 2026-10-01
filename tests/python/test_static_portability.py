from pathlib import Path


def test_active_code_has_no_machine_specific_paths():
    root = Path(__file__).resolve().parents[2]
    forbidden = ("/Users/", "/sci/", "MY_R", "Zaslab/entropy", "Zaslab/Entropy")
    files = []
    for directory in ("src", "scripts", "slurm", "explorer"):
        files.extend(path for path in (root / directory).rglob("*") if path.is_file() and "__pycache__" not in path.parts)
    failures = {
        str(path.relative_to(root)): token
        for path in files
        for token in forbidden
        if token in path.read_text(errors="ignore")
    }
    assert not failures

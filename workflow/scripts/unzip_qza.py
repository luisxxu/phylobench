import argparse
import zipfile
from pathlib import Path


def extract_members(qza: Path, out_dir: Path, names: list[str]) -> None:
    with zipfile.ZipFile(qza) as archive:
        roots = {n.split("/", 1)[0] for n in archive.namelist()}
        if len(roots) != 1:
            raise ValueError(f"{qza}: expected one archive root, found {sorted(roots)}")
        (root,) = roots
        out_dir.mkdir(parents=True, exist_ok=True)
        for name in names:
            target = (out_dir / name).resolve()
            if target.parent != out_dir.resolve():
                raise ValueError(f"member name escapes the output dir: {name}")
            target.write_bytes(archive.read(f"{root}/data/{name}"))


def main():
    parser = argparse.ArgumentParser(description="Extract data/ members of a .qza.")
    parser.add_argument("qza", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("members", nargs="+")
    args = parser.parse_args()
    extract_members(args.qza, args.out_dir, args.members)


if __name__ == "__main__":
    main()

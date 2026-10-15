import argparse
from pathlib import Path

from phylobench.simulate import tiny_fixture


def main():
    parser = argparse.ArgumentParser(description="Regenerate the tiny test fixture.")
    parser.add_argument("out_dir", type=Path)
    out_dir = parser.parse_args().out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, text in tiny_fixture().items():
        (out_dir / name).write_text(text)


if __name__ == "__main__":
    main()

import json
import sys
from pathlib import Path

from phylobench.jplace import from_sepp


def main():
    source, out = (Path(arg) for arg in sys.argv[1:3])
    out.write_text(json.dumps(from_sepp(json.loads(source.read_text()))))


if __name__ == "__main__":
    main()

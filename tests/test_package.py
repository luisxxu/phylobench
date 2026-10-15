import importlib


def test_package_imports():
    assert importlib.import_module("phylobench").__name__ == "phylobench"

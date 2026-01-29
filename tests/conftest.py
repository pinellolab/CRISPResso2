import contextlib
import os
import sys
import types
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "mplconfig"))


try:
    import pytest_check  # noqa: F401
except Exception:
    fallback = types.ModuleType("pytest_check")

    class _CheckContext:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def equal(self, a, b):
            assert a == b

        def is_true(self, value):
            assert value

        def is_false(self, value):
            assert not value

    fallback.check = _CheckContext()
    sys.modules["pytest_check"] = fallback


try:
    from CRISPResso2 import CRISPResso2Align  # noqa: F401
except Exception as exc:
    raise pytest.UsageError(
        "CRISPResso2 C-extensions are not importable in this Python. "
        "Rebuild with the same interpreter, e.g.:\n"
        "  python setup.py build_ext --inplace\n"
        f"Original import error: {exc}"
    )

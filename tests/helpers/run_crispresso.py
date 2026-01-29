import importlib
import os
import sys


def _patch_sysconf():
    original_sysconf = os.sysconf

    def safe_sysconf(name):
        if name == 'SC_SEM_NSEMS_MAX':
            return 256
        return original_sysconf(name)

    os.sysconf = safe_sysconf


def main():
    if len(sys.argv) < 2:
        raise SystemExit("Usage: run_crispresso.py <ModuleName> -- <args...>")

    if '--' in sys.argv:
        split_idx = sys.argv.index('--')
        module_name = sys.argv[1]
        args = sys.argv[split_idx + 1:]
    else:
        module_name = sys.argv[1]
        args = sys.argv[2:]

    _patch_sysconf()

    module = importlib.import_module(f"CRISPResso2.{module_name}")
    sys.argv = [module_name] + args
    module.main()


if __name__ == '__main__':
    main()

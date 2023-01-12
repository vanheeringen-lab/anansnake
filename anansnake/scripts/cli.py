#!/usr/bin/env python

import subprocess as sp
from os.path import dirname, join
import sys
from anansnake import __version__
from seq2science import __version__ as version_s2s  # noqa
from snakemake import __version__ as version_sm


def cli():
    """
    wrap snakemake with the anansnake Snakefile preset
    """
    script_dir = dirname(__file__)
    pkg_dir = dirname(script_dir)
    snakefile = join(pkg_dir, "Snakefile")
    cmd = " ".join([f"snakemake --snakefile {snakefile} --use-conda --conda-frontend mamba"] + sys.argv[1:])

    if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
        print(f"anansnake v{__version__}")
        print(f"seq2science v{version_s2s}")
        print(f"snakemake v{version_sm}")
        return

    retcode = sp.call(cmd, shell=True)
    print("")  # newline
    if retcode != 0:
        sys.exit(retcode)


if __name__ == "__main__":
    cli()

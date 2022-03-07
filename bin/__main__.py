
# Imports
import argparse


def minetwister():
    parser = argparse.ArgumentParser(
        description="Mine genomes for twister ribozymes"
    )

    parser.add_argument(
        "-i", dest="reference", help="path to reference genome", required=True
    )

    args = parser.parse_args()
    print(args)
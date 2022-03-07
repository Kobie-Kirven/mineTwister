
# Imports
import argparse
import subprocess
from minetwister import *

def minetwister():
    parser = argparse.ArgumentParser(
        description="Mine genomes for twister ribozymes"
    )

    parser.add_argument(
        "-i", dest="reference", help="path to reference genome", required=True
    )

    parser.add_argument(
        "-c", dest="query", help="concensus twister sequence", required=True
    )

    parser.add_argument(
        "-o", dest="output", help="output", required=True
    )

    args = parser.parse_args()
    print(args)

    # Run minetwister
    build_blast_db(args.reference)
    blast_query(args.query, args.reference, args.output)

if __name__ == "__main__":
    minetwister()
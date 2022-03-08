################################################################################
# Mine a genome of interest for twister ribozymes
#
# Author: Kobie Kirven
################################################################################


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
        "-s", dest="singularity", help="singularity path", required=True
    )

    parser.add_argument(
        "-d", dest="data", help="singularity data path", required=True
    )

    parser.add_argument(
        "-o", dest="output", help="output", required=True
    )



    args = parser.parse_args()
    print(args)

    #---------------------------------------------------------------------------
    #  Run minetwister
    #---------------------------------------------------------------------------
    
    # Build blast database from reference genome
    build_blast_db(args.reference)

    # Blast query against reference genome
    blast_query(args.query, args.reference, "blast_output.tab")
    
    # Parse blast output
    hits = parse_blast_output("blast_output.tab")
    with open("blast_fasta.fasta", "w") as fh:
        for i, hit in enumerate(hits):
            if hits[-1]== "Minus":
                fh.write(f">potential_twister_{i}"+"\n"+rev_transcribe(get_flanking_seq(args.reference,hit[8], hit[9],50))+
                "\n")
            else:
                fh.write(f">potential_twister_{i}" + "\n" + get_flanking_seq(args.reference,hit[8], hit[9],50) + "\n")

    # Run r2dt
    run_r2dt(args.data, args.singularity)

if __name__ == "__main__":
    minetwister()
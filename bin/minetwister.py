###############################################################################
# Functions for minetwister
###############################################################################
import subprocess

def build_blast_db(reference):
    """
    Build blast database from reference genome
    """
    print("Building blast database from reference genome")
    subprocess.run(
        ["makeblastdb", "-in", reference, "-dbtype", "nucl"]
    )


def blast_query(query, reference, output):
    """
    Blast query against reference genome
    """
    print("Blasting query against reference genome")
    subprocess.run(
        ["blastn", "-query", query, "-db", reference, "-out", output]
    )


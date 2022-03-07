###############################################################################
# Functions for minetwister
###############################################################################
import subprocess
from Bio import SeqIO

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
        ["blastn", "-query", query, "-db", reference, "-out", output, "-outfmt", "6"]
    )


def parse_blast_output(blast_output):
    """
    Parse blast output
    """
    return [line.strip("\n").split("\t") for line in open(blast_output).readlines()]



def get_flanking_seq(genome, start, end, flanking_length):
    """
    Get flanking based on Blast hit
    """
    for rec in SeqIO.parse(genome, "fasta"):
        return str(rec.seq[start-flanking_length:end+flanking_length])
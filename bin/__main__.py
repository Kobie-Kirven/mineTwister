################################################################################
# Mine a genome of interest for twister ribozymes
#
# Author: Kobie Kirven
################################################################################


# Imports
import argparse
import subprocess
from Bio import SeqIO

###############################################################################
# Functions for minetwister
###############################################################################

def build_blast_db(reference):
    """
    Build blast database from reference genome
    """
    print("Building blast database from reference genome")
    p1 = subprocess.run(
        ["makeblastdb", "-in", reference, "-dbtype", "nucl"]
    )
    p1.wait()


def blast_query(query, reference, output):
    """
    Blast query against reference genome
    """
    print("Blasting query against reference genome")

    extended_command = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand"
    p1 = subprocess.run(
        ["blastn", "-query", query, "-db", reference, "-out", output, "-outfmt",extended_command ]
    )
    p1.wait()

def parse_blast_output(blast_output):
    """
    Parse blast output
    """
    return [line.strip("\n").split("\t") for line in open(blast_output).readlines()]

def rev_transcribe(dna):
    """
    Transcribe DNA to RNA
    """
    replace_dict = {"T":"A", "G":"C", "A":"T", "C":"G"}
    return "".join([replace_dict[base] for base in dna][::-1])


def get_flanking_seq(genome, start, end, flanking_length):
    """
    Get flanking based on Blast hit
    """
    for rec in SeqIO.parse(genome, "fasta"):
        return str(rec.seq[start-flanking_length:end+flanking_length])


def run_r2dt(cms_data_path, singularity_image_path):
    p1 = subprocess.run(["singularity","exec", "--bind",
    cms_data_path + ":/rna/r2dt/data/cms",singularity_image_path, "r2dt.py", "draw", 
    "blast_fasta.fasta", "./temp_res"])
    p1.wait()

def build_html_output(genome, scaffold, start, end, sequence, figure_path):
    """
    Build html output
    """
    html_output = """
    <html>
    <head>
    <title>MineTwister</title>
    </head>
    <body>
    <h1>Potential Twisters</h1>
    <p>
    <h3>Genome: {}</h3>
    <h3>Scaffold: {}</h3>
    <h3>Start: {}</h3>
    <h3>End: {}</h3>
    <h3>Sequence: {}</h3>
    </p>
    <img src="{}" alt="image">
    </body>
    </html>
    """.format(genome, scaffold, start, end, sequence, figure_path)
    return html_output

###############################################################################
# mineTwister main
###############################################################################
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
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

    extended_command = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand"
    subprocess.run(
        ["blastn", "-query", query, "-db", reference, "-out", output, "-outfmt",extended_command ]
    )


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
    subprocess.run(["singularity","exec", "--bind",
    cms_data_path + ":/rna/r2dt/data/cms",singularity_image_path, "r2dt.py", "draw", 
    "blast_fasta.fasta", "./temp_res"])


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
import subprocess
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO 
from io import StringIO
from Bio import SeqIO

def RunDiamond(query, db, output, qcov = 80, maxseq = 5, threads=12):
    """Runs a DIAMOND search.
    Args:
        query (str): Path to the query FASTA file.
        db (str): Path to the DIAMOND database.
        output (str): Path to the output file.
        qcov (int, optional): Query coverage threshold. Defaults to 80.
        maxseq (int, optional): Maximum number of target sequences to report. Defaults to 5.
        threads (int, optional): Number of threads to use. Defaults to 12.
    """
    command = f"diamond blastp -d {db} -q {query} -o {output} -p {threads} --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore ppos full_qseq full_sseq -b4 --query-cover {qcov} -k {maxseq} --no-self-hits"
    subprocess.run(command, shell=True)

def ProteinAligner(fasta):
    """Aligns a FASTA file of protein sequences using MAFFT.
    Args:        fasta (str): Path to the input FASTA file.
    Returns:        AlignIO.MultipleSeqAlignment: The aligned sequences.
    """
    fasta = "\n".join(fasta)
    mafft_cline = MafftCommandline(input="-")
    stdout, stderr = mafft_cline(stdin=fasta)
    ProteinAlignment = AlignIO.read(StringIO(stdout), "fasta")
    return ProteinAlignment

def RunCDHIT(input, output, max_identity=0.95, min_cov = 0.8):
    """Runs CD-HIT to cluster sequences in a FASTA file based on sequence identity.
    Args:        input (str): Path to the input FASTA file.
        output (str): Path to the output FASTA file where the clustered sequences will be saved.
        max_identity (float, optional): The maximum sequence identity threshold for clustering. Defaults to 0.95.
        min_cov (float, optional): The minimum coverage threshold for clustering. Defaults to 0.8."""
    command = [
        "cd-hit",
        "-i", input,
        "-o", output,
        "-c", str(max_identity),
        "-T", str(12),
        "-p", str(1),
        "-sc", str(1),
        "-sf", str(1),
        "-aL" ,str(min_cov),
    ]
    subprocess.run(command, check=True)
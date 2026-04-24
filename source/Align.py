import subprocess
def RunDiamond(query, db, output, qcov = 80, maxseq = 5, threads=4):
    """Runs a DIAMOND search.
    Args:
        query (str): Path to the query FASTA file.
        db (str): Path to the DIAMOND database.
        output (str): Path to the output file.
        qcov (int, optional): Query coverage threshold. Defaults to 80.
        maxseq (int, optional): Maximum number of target sequences to report. Defaults to 5.
        threads (int, optional): Number of threads to use. Defaults to 4.
    """
    command = f"diamond blastp -d {db} -q {query} -o {output} -p {threads} --outfmt 6 -b4 --query-cover {qcov} -k {maxseq}"
    subprocess.run(command, shell=True)
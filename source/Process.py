import sys
sys.path.append('./')
import subprocess
from Bio import SeqIO
from source import Align

def InsertTag(input, output, tag):
    """Insert a database-refering tag at the end of each sequence header in a FASTA file.
    Args:
        input: the input FASTA file.
        output: the output FASTA file.
        tag: the tag to be inserted at the end of each sequence header.
    Returns:
        Creates a new FASTA file with the tag inserted at the end of each sequence header.
        Creastes a new file with the list of all sequence headers in the input FASTA file, without the tag.
    """
    Counter = 0
    with open(input, 'r') as f:
        lines = f.readlines()
    with open(output, 'w+') as f:
        for line in lines:
            if line.startswith('>'):
                line = line.lstrip('>')
                f.write(f">{tag}_{Counter}|{tag}\n")
                Counter += 1
            else:
                f.write(line)

    Counter = 0
    with open(f"{output}_ids", 'w+') as f:
        for line in lines:
            if line.startswith('>'):
                line = line.lstrip('>')
                f.write(f">{tag}_{Counter}|{line}\n")
                Counter += 1
            else:
                pass
def ConcatenateFiles(path,output,extension = "fasta"):
    """Concatenate all fasta files in a given path into a single output file.
    Args:
        path (str): Path to the directory containing the files to concatenate.
        output (str): Path to the output file where the concatenated content will be saved.
        extension (str): The file extension to look for (default is "fasta").
    """
    command = f"cat {path}/*.{extension} > {output}"
    subprocess.run(command, shell=True, check=True)
def CreateDatabase(input, output):
    """Creates a DIAMOND database from a FASTA file.
    Args:
        input (str): Path to the input FASTA file.
        output (str): Path to the output DIAMOND database.
    """
    command = f"diamond makedb --in {input} -d {output}"
    subprocess.run(command, shell=True)
def GetSequences(fasta):
    """Reads a FASTA file and returns a list of sequences.
    Args:        fasta (str): Path to the input FASTA file.
    Returns:        list: A list of sequences.
    """
    with open(fasta, "r") as handle:
        sequences = [f">{record.id}\n{record.seq}" for record in SeqIO.parse(handle, "fasta")]
    return sequences
def CountId(path):
    """Counts the number of sequences in a FASTA file.
    Args:
        path (str): Path to the FASTA file.
    Returns:
        int: The number of sequences in the FASTA file.
    """
    Command = ["grep", "-c", ">", path]
    result = subprocess.run(Command, capture_output=True, text=True, check=True)
    return int(result.stdout.strip())
def ClusterCounter(path):
    """Runs CD-HIT at different identity thresholds and counts the number of clusters formed at each threshold.
    Args:
        path (str): Path to the input FASTA file.
    Returns:
        tuple: A tuple containing two lists: the first list contains the number of clusters at each threshold, and the second list contains the corresponding percentages of the total number of sequences.
    """
    TotalSize = CountId(path)
    counts = []
    perc = []
    for i in [1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7]:
        Align.RunCDHIT(
            input = path,
            output = f"{path}.out",
            max_identity = i
        )
        PartialSize = CountId(f"{path}.out")
        Percentage = (PartialSize / TotalSize) * 100
        counts.append(PartialSize)
        perc.append(Percentage)
        subprocess.run(["rm", f"{path}.out"], check=True)
        subprocess.run(["rm", f"{path}.out.clstr"], check=True)
    return counts, perc
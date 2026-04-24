import subprocess

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

def ConcatenateFiles(path,output):
    """Concatenate all fasta files in a given path into a single output file.
    Args:
        path (str): Path to the directory containing the files to concatenate.
        output (str): Path to the output file where the concatenated content will be saved.
    """
    command = f"cat {path}/*.fasta > {output}"
    subprocess.run(command, shell=True, check=True)

def CreateDatabase(input, output):
    """Creates a DIAMOND database from a FASTA file.
    Args:
        input (str): Path to the input FASTA file.
        output (str): Path to the output DIAMOND database.
    """
    command = f"diamond makedb --in {input} -d {output}"
    subprocess.run(command, shell=True)
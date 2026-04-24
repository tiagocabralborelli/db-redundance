import subprocess

def ToProtein(input,output,minlength):
    """Converts a fasta file to a protein fasta file using the translate function from the source module
    Args:
        input (str): The path to the input fasta file
        output (str): The path to the output fasta file
        minlength (int): The minimum length of the protein sequences to be included in the output file
    Requires:
        seqkit toolkit
    """
    command = f"seqkit translate {input} -T 11 -s -m {minlength} > {output}"
    subprocess.run(command, shell=True, check=True)
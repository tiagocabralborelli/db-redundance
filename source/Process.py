import sys
sys.path.append('./')
import subprocess
from Bio import SeqIO
from source import Align
import json

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

def CreateMetadataFile(Input, SCHEMAS):
    """Creates a metadata file in JSON format from a FASTA file containing sequence headers with database-specific metadata.
    Args:
        Input (str): Path to the input FASTA file containing sequence headers with database-specific metadata.
        SCHEMAS (dict): A dictionary containing the split points for each database's metadata schema and the corresponding index information.
    Returns:        dict: A dictionary containing the converged metadata, with keys "Drug Class" and "Name".
    """
    MetadataDict = {}
    with open(Input, "r") as file:
        for line in file:
            if line.startswith(">"):
                ID = line.strip(">").strip().split("|")[0]
                ## Check database specificity
                DatabaseTag = ID.split("_")[0]
                MetadataDict[ID] = ConvergeSchemas(DatabaseTag, line, SCHEMAS)

    return MetadataDict

def ConvergeSchemas(Database, String, SCHEMAS):
    """Converges the different metadata schemas of the databases into a single, unified schema.
    Args:        
        Database (str): The name of the database from which the metadata is being extracted.
        String (str): The string containing the metadata to be extracted.
        SCHEMAS (dict): A dictionary containing the split points for each database's metadata schema and the corresponding index information.
    Returns:        dict: A dictionary containing the converged metadata, with keys "Drug Class" and "Name".
    """

    match Database:
        case "CARD":
            ARO = String.strip(">").strip().split("|")[SCHEMAS["CARD"]["AROSplitPoint"]]
            DrugClass = SCHEMAS["CARD"]["IndexInfo"][ARO]["Drug Class"]
            Name = SCHEMAS["CARD"]["IndexInfo"][ARO]["ARO Name"]
            return {"Drug Class": DrugClass, "Name": Name}

        case "NDARO":
            RefSeq = String.strip(">").strip().split("|")[-1].split(" ")[SCHEMAS["NDARO"]["AccSplitPoint"]]
            try:
                DrugClass = SCHEMAS["NDARO"]["IndexInfo"][RefSeq]["Class"]
            except:
                DrugClass = "Not Found"
            try:
                Name = SCHEMAS["NDARO"]["IndexInfo"][RefSeq]["Gene family"]
            except:
                Name = "Not Found"
            return {"Drug Class": DrugClass.strip().lower(), "Name": Name}

        case "MEGARES":
            DrugClass = String.strip(">").strip().split("|")[SCHEMAS["MEGARES"]["DrugClassSplitPoint"]]
            Name  = String.strip(">").strip().split("|")[SCHEMAS["MEGARES"]["NameSplitPoint"]]
            return {"Drug Class": DrugClass, "Name": Name.split("_")[0]}
        
        case "HMD":
            DrugClass = String.strip(">").strip().split("|")[SCHEMAS["HMD"]["DrugClassSplitPoint"]]
            Name  = String.strip(">").strip().split("|")[SCHEMAS["HMD"]["NameSplitPoint"]]
            return {"Drug Class": DrugClass, "Name": Name}
        
        case "NCRD":
            DrugClass = String.strip(">").strip().split("|")[SCHEMAS["NCRD"]["DrugClassSplitPoint"]]
            Name = String.strip(">").strip().split("|")[SCHEMAS["NCRD"]["NameSplitPoint"]]
            return {"Drug Class": DrugClass, "Name": Name}
        
        case "RESFINDER":            
            DrugClass = String.strip(">").strip().split("|")[SCHEMAS["RESFINDER"]["DrugClassSplitPoint"]]
            Name = String.strip(">").strip().split("|")[SCHEMAS["RESFINDER"]["NameSplitPoint"]]
            return {"Drug Class": DrugClass.split("_")[0], "Name": Name.split("_")[0]}
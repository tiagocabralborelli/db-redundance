# An analysis of the antimicrobial resistance reference database redundancy and its impact on functional annotation

## Background

## Documentation 
### Modules and fuctions
#### Translate
```
ToProtein(input,output,minlength):
    Converts a fasta file to a protein fasta file using the translate function from the source module
    Args:
        input (str): The path to the input fasta file
        output (str): The path to the output fasta file
        minlength (int): The minimum length of the protein sequences to be included in the output file
    Requires:
        seqkit toolkit
```
#### Process 
```
InsertTag(input, output, tag):
    Insert a database-refering tag at the end of each sequence header in a FASTA file.
    Args:
        input: the input FASTA file.
        output: the output FASTA file.
        tag: the tag to be inserted at the end of each sequence header.
    Returns:
        Creates a new FASTA file with the tag inserted at the end of each sequence header.
        Creates a new file with the list of all sequence headers in the input FASTA file, without the tag.

ConcatenateFiles(path,output,extension):
    Concatenate all fasta files in a given path into a single output file.
    Args:
        path (str): Path to the directory containing the files to concatenate.
        output (str): Path to the output file where the concatenated content will be saved.
        extension (str): The file extension to look for (default is "fasta").

CreateDatabase(input, output):
    Creates a DIAMOND database from a FASTA file.
    Args:
        input (str): Path to the input FASTA file.

GetSequences(fasta):
    Reads a FASTA file and returns a list of sequences.
    Args:        
        fasta (str): Path to the input FASTA file.
    Returns:       
        list: A list of sequences.

ClusterCounter(path):
    Runs CD-HIT at different identity thresholds and counts the number of clusters formed at each threshold.
    Args:
        path (str): Path to the input FASTA file.
    Returns:
        tuple: A tuple containing two lists: the first list contains the number of clusters at each threshold, and the second list contains the corresponding percentages of the total number of sequences.

def CreateMetadataFile(Input, Output, SCHEMAS):
    Creates a metadata file in JSON format from a FASTA file containing sequence headers with database-specific metadata.
    Args:
        Input (str): Path to the input FASTA file containing sequence headers with database-specific metadata.
        Output (str): Path to the output JSON file where the metadata will be saved.
        SCHEMAS (dict): A dictionary containing the split points for each database's metadata schema and the corresponding index information.
    Returns:        dict: A dictionary containing the converged metadata, with keys "Drug Class" and "Name".
```
#### Align
```
RunDiamond(query, db, output, qcov = 80, maxseq = 5, threads=4):
    Runs a DIAMOND search.
    Args:
        query (str): Path to the query FASTA file.
        db (str): Path to the DIAMOND database.
        output (str): Path to the output file.
        qcov (int, optional): Query coverage threshold. Defaults to 80.
        maxseq (int, optional): Maximum number of target sequences to report. Defaults to 5.
        threads (int, optional): Number of threads to use. Defaults to 4.
ProteinAligner(fasta):
    Aligns a FASTA file of protein sequences using MAFFT.
    Args:        
        fasta (str): Path to the input FASTA file.
    Returns:        
        AlignIO.MultipleSeqAlignment: The aligned sequences.
RunCDHIT(input, output, max_identity=0.95, min_cov = 0.8):
    Runs CD-HIT to cluster sequences in a FASTA file based on sequence identity.
    Args:        input (str): Path to the input FASTA file.
        output (str): Path to the output FASTA file where the clustered sequences will be saved.
        max_identity (float, optional): The maximum sequence identity threshold for clustering. Defaults to 0.95.
        min_cov (float, optional): The minimum coverage threshold for clustering. Defaults to 0.8."""
```
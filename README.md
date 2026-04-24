# An analysis of the antimicrobial resistance reference database redundancy and its impact on functional annotation



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
#### Standardize 
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

ConcatenateFiles(path,output):
    Concatenate all fasta files in a given path into a single output file.
    Args:
        path (str): Path to the directory containing the files to concatenate.
        output (str): Path to the output file where the concatenated content will be saved.

CreateDatabase(input, output):
    Creates a DIAMOND database from a FASTA file.
    Args:
        input (str): Path to the input FASTA file.
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
```
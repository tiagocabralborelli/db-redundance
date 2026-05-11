from pathlib import Path
import subprocess

def MakeDiamondDB (fasta_dir, db_dir):
    """
    Concatenate all FASTA files from a folder and create a single DIAMOND database.
    
    Args:
        fasta_dir (str or Path): Folder containing input FASTA files.
        db_dir (str or Path): Output folder for the concatenated FASTA and DIAMOND database.
    """
    fasta_dir = Path(fasta_dir)
    db_dir = Path(db_dir)
    
    db_dir.mkdir(parents=True, exist_ok=True)

    output_fasta = db_dir / "combined.fasta"
    output_db = db_dir / "combined"

    fasta_files = sorted(fasta_dir.glob("*.fasta"))

    if len(fasta_files) == 0:
        raise FileNotFoundError(f"No FASTA files found in {fasta_dir}")
    
    with open(output_fasta, "w") as outfile:
        for fasta_file in fasta_files:
            print(f"Adding: {fasta_file.name}")
            with open(fasta_file, "r") as infile:
                outfile.write(infile.read())
            outfile.write("\n")  # Ensure separation between files

    print(f"Concatenated FASTA saved to: {output_fasta}")

    command = [
        "diamond", "makedb",
        "--in", str(output_fasta),
        "--db", str(output_db)
    ]

    print ("Running:", " ".join(command))

    with subprocess.Popen(command, stderr=subprocess.PIPE, text=True) as process:
        for line in process.stderr:
            print(line, end="")

    if process.returncode != 0:
        raise RuntimeError(f"DIAMOND makedb failed with return code {process.returncode}")
    
    print (f"DIAMOND database saved as: {output_db}.dmnd")


def PrepareGenomeQueryFasta(genomes_dir, output_filename="all_genomes_query.fasta"):
    """
    Rename headers of protein FASTA files from NCBI genome folders and concatenate
    into a single query FASTA file.

    Args:
        genomes_dir (str or Path): Folder containing one subfolder per genome (GCA_...).
        output_filename (str): Name of the output concatenated FASTA file.
    """
    genomes_dir = Path(genomes_dir)
    output_fasta = genomes_dir / output_filename

    genome_folders = sorted([f for f in genomes_dir.iterdir() if f.is_dir()])

    if len(genome_folders) == 0:
        raise FileNotFoundError(f"No genome folders found in {genomes_dir}")

    with open(output_fasta, "w") as outfile:
        for genome_folder in genome_folders:
            gca_id = genome_folder.name
            protein_file = genome_folder / "protein.faa"

            if not protein_file.exists():
                print(f"Warning: no protein.faa found in {gca_id}, skipping.")
                continue

            print(f"Processing: {gca_id}")

            with open(protein_file, "r") as infile:
                for line in infile:
                    if line.startswith(">"):
                        protein_id = line[1:].split()[0]
                        rest = " ".join(line[1:].split()[1:])
                        outfile.write(f">{gca_id}|{protein_id} {rest}\n")
                    else:
                        outfile.write(line)

            outfile.write("\n")

    print(f"Query FASTA saved to: {output_fasta}")
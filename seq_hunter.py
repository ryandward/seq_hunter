import datetime
import sys
import gzip
import zlib
import argparse
import configparser
import time
import os
from urllib.request import urlretrieve
from urllib.error import HTTPError, URLError
from Bio import Entrez
from Bio import SeqIO
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.console import Console

DESCRIPTION_WIDTH = 40
CONFIG_FILE = "config.ini"
console = Console()

def get_email():
    """Retrieve or prompt for the user's email for NCBI Entrez API."""
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)

    if not config.has_section("Entrez"):
        config.add_section("Entrez")

    if not config.has_option("Entrez", "email"):
        console.print("NCBI requires an email address for its API usage.", style="cyan")
        email = input("Please enter your email address: ")
        config.set("Entrez", "email", email)

        with open(CONFIG_FILE, "w") as configfile:
            config.write(configfile)
    else:
        email = config.get("Entrez", "email")
        console.print(f"Using stored email address: {email}", style="cyan")

    return email

# Set your email for NCBI Entrez
Entrez.email = get_email()

def download_sequences_by_accession(accession_numbers, source):
    """Download sequences from NCBI by accession numbers and the specified source."""
    assembly_data = fetch_assembly_data(accession_numbers, source)
    assemblies_sequences = download_assemblies_sequences(assembly_data, source)
    return assemblies_sequences

def fetch_assembly_data(accession_numbers, source):
    """Fetch assembly data for the given accession numbers and source."""
    assembly_data = []

    with Progress(
        TextColumn(
            f"[bold cyan]{{task.description:<{DESCRIPTION_WIDTH}}}", justify="left"
        ),
        BarColumn(bar_width=None),
        TextColumn("{task.percentage:>3.0f}%", justify="right"),
        transient=False,
    ) as progress:
        fetch_task = progress.add_task(
            "[1/3] Fetching assembly accessions...", total=len(accession_numbers)
        )

        for accession_number in accession_numbers:
            db = "nuccore" if len(accession_number) == 11 else "assembly"
            handle = Entrez.esearch(db=db, term=accession_number)
            record = Entrez.read(handle)

            if not record["IdList"]:
                console.print(
                    f"[bold red]Accession number not found: {accession_number}"
                )
                progress.update(fetch_task, advance=1, refresh=True)
                continue

           
            if db == "nuccore":
                assembly_id = record["IdList"][0]
                try:
                    handle = Entrez.esummary(db="assembly", id=assembly_id)
                except (HTTPError, URLError) as e:
                    console.print(
                        f"[bold red]Error while fetching assembly summary for {assembly_id}: {str(e)}"
                    )
                    continue

                try:
                    summary = Entrez.read(handle)
                except RuntimeError as e:
                    console.print(
                        f"[bold red]Error while reading assembly summary for {assembly_id}: {str(e)}"
                    )
                    continue

                assembly_data.append({
                    "accession": summary["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"],
                    "ftp_path": summary["DocumentSummarySet"]["DocumentSummary"][0][f"FtpPath_{source}"]
                })
            else:
                for assembly_id in record["IdList"]:
                    try:
                        handle = Entrez.esummary(db="assembly", id=assembly_id)
                    except (HTTPError, URLError) as e:
                        console.print(
                            f"[bold red]Error while fetching assembly summary for {assembly_id}: {str(e)}"
                        )
                        continue

                    try:
                        summary = Entrez.read(handle)
                    except RuntimeError as e:
                        console.print(
                            f"[bold red]Error while reading assembly summary for {assembly_id}: {str(e)}"
                        )
                        continue

                    assembly_data.append({
                        "accession": summary["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"],
                        "ftp_path": summary["DocumentSummarySet"]["DocumentSummary"][0][f"FtpPath_{source}"]
                    })

            console.print(summary["DocumentSummarySet"]["DocumentSummary"][0]["Synonym"]["Genbank"])
            progress.update(fetch_task, advance=1, refresh=True)

    return assembly_data

def download_assemblies_sequences(assembly_data, source):
    """Download sequences for the given assembly data."""
    assemblies_sequences = {}

    with Progress(
        TextColumn(
            f"[bold cyan]{{task.description:<{DESCRIPTION_WIDTH}}}", justify="left"
        ),
        BarColumn(bar_width=None),
        TextColumn("{task.percentage:>3.0f}%", justify="right"),
        transient=False,
    ) as progress:
        download_task = progress.add_task(
            "[2/3] Downloading sequences...", total=len(assembly_data)
        )

        for assembly in assembly_data:
            ftp_path = assembly["ftp_path"]
            if ftp_path:
                genbank_file_name = ftp_path.split("/")[-1] + f"_genomic.gbff.gz"
                genbank_url = f"{ftp_path}/{genbank_file_name}"

                try:
                    compressed_file, _ = urlretrieve(genbank_url)
                except URLError as e:
                    console.print(
                        f"[bold red]Error while downloading sequence for accession {assembly['accession']}: {str(e)}"
                    )
                    continue

                try:
                    with gzip.open(compressed_file, 'rt') as compressed_handle:
                        try:
                            sequences = list(SeqIO.parse(compressed_handle, "genbank"))
                        except ValueError as e:
                            console.print(
                                f"[bold red]Error while parsing sequence for accession {assembly['accession']}: {str(e)}"
                            )
                            continue

                        assemblies_sequences[assembly['accession']] = sequences
                except zlib.error as e:
                    console.print(
                        f"[bold red]Error while opening gzipped file for accession {assembly['accession']}: {str(e)}"
                    )
                    continue

            progress.update(download_task, advance=1, refresh=True)

    return assemblies_sequences


def save_sequences_to_file(assemblies_sequences):
    """Save the downloaded sequences to individual files."""
    total_sequences = sum(len(sequences) for sequences in assemblies_sequences.values())

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    assemblies_file_name = f"assemblies_{timestamp}.tsv"

    with Progress(
        TextColumn(
            f"[bold cyan]{{task.description:<{DESCRIPTION_WIDTH}}}", justify="left"
        ),
        BarColumn(bar_width=None),
        TextColumn("{task.percentage:>3.0f}%", justify="right"),
        transient=False,
    ) as progress:
        save_task = progress.add_task(
            "[3/3] Saving sequences to files...", total=len(assemblies_sequences)
        )

        # Write header to organisms file
        with open(assemblies_file_name, "w") as assemblies_file:
            assemblies_file.write("AssemblyAccession\tOrganism\n")

        for accession, sequences in assemblies_sequences.items():
            os.makedirs("sequences", exist_ok=True)

            file_path = f"sequences/{accession}.gb"

            try:
                with open(file_path, "w") as output_handle:
                    SeqIO.write(sequences, output_handle, "genbank")
            except IOError as e:
                console.print(
                    f"[bold red]Error while writing sequence to file {file_path}: {str(e)}"
                )
                continue

            progress.update(save_task, advance=1, refresh=True)

            # Check that all organisms in the assembly are the same
            organisms = {seq.annotations["organism"] for seq in sequences}
            if len(organisms) == 1:
                organism = organisms.pop()
                with open(assemblies_file_name, "a") as assemblies_file:
                    assemblies_file.write(f"{accession}\t{organism}\n")
    
    return assemblies_file_name


def main():
    parser = argparse.ArgumentParser(
        description="Download sequences from NCBI by accession numbers."
    )
    parser.add_argument(
        "accession_numbers",
        metavar="ACCESSION",
        nargs="+",
        help="An accession number of a sequence to download.",
    )
    parser.add_argument(
        "-s",
        "--source",
        choices=["GenBank", "RefSeq"],
        default="GenBank",
        help="The source database of the sequences (default: %(default)s).",
    )

    args = parser.parse_args()

    try:
        assemblies_sequences = download_sequences_by_accession(
            args.accession_numbers, args.source
        )
        assemblies_file_name = save_sequences_to_file(assemblies_sequences)
    except KeyboardInterrupt:
        console.print("\n[bold red]Interrupted by user. Exiting.")
        sys.exit(1)

    console.print(
        "[bold green]Script completed successfully.",
        style="bold green",
    )
    console.print(
        "Please check the 'sequences' folder for the downloaded sequences,",
        f"and '{assemblies_file_name}' for a list of the accession numbers and organisms.",
        style="bold yellow",
    )

if __name__ == "__main__":
    main()

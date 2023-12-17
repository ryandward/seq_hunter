import datetime
import sys
import gzip
import zlib
import argparse
import configparser
import os
from urllib.request import urlretrieve
from urllib.error import HTTPError, URLError
from Bio import Entrez
from Bio import SeqIO
from rich.progress import Progress, BarColumn, TextColumn
from rich.console import Console
from rich.table import Table

DESCRIPTION_WIDTH = 40
CONFIG_FILE = "config.ini"
DOWNLOAD_LOG_FILE = "download_log.tsv"
console = Console()


def get_email():
    """Retrieve or prompt for the user's email for NCBI Entrez API."""
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)

    if not config.has_section("Entrez"):
        config.add_section("Entrez")

    email = config.get("Entrez", "email", fallback=None)
    if email is None:
        console.print(
            "NCBI requires an email address for its API usage.", style="cyan")
        email = input("Please enter your email address: ")
        config.set("Entrez", "email", email)
        with open(CONFIG_FILE, "w") as configfile:
            config.write(configfile)
    else:
        console.print(f"Using stored email address: {email}", style="cyan")

    return email


def fetch_assembly_data(accession_numbers, source, progress):
    """Fetch assembly data for the given accession numbers and source."""
    assembly_data = []
    fetch_task = progress.add_task(
        "[1/3] Fetching assembly accessions...", total=len(accession_numbers))

    for accession_number in accession_numbers:
        db = "nuccore" if len(accession_number) == 11 else "assembly"
        handle = Entrez.esearch(db=db, term=accession_number)
        record = Entrez.read(handle)

        if not record["IdList"]:
            console.print(
                f"[bold red]Accession number not found: {accession_number}")
            progress.update(fetch_task, advance=1, refresh=True)
            continue

        assembly_id = record["IdList"][0]
        try:
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(handle)
        except (HTTPError, URLError, RuntimeError) as e:
            console.print(
                f"[bold red]Error while processing {assembly_id}: {str(e)}")
            continue

        accession_key = "Genbank" if source == "GenBank" else "RefSeq"
        assembly_data.append({
            "accession": summary["DocumentSummarySet"]["DocumentSummary"][0]["Synonym"][accession_key],
            "ftp_path": summary["DocumentSummarySet"]["DocumentSummary"][0][f"FtpPath_{source}"]
        })
        progress.update(fetch_task, advance=1, refresh=True)

    return assembly_data


def download_assemblies_sequences(assembly_data, source, progress):
    """Download sequences for the given assembly data."""
    assemblies_sequences = {}
    download_task = progress.add_task(
        "[2/3] Downloading sequences...", total=len(assembly_data))

    for assembly in assembly_data:
        ftp_path = assembly["ftp_path"]
        if not ftp_path:
            continue

        genbank_file_name = ftp_path.split("/")[-1] + f"_genomic.gbff.gz"
        genbank_url = f"{ftp_path}/{genbank_file_name}"

        try:
            compressed_file, _ = urlretrieve(genbank_url)
            with gzip.open(compressed_file, 'rt') as compressed_handle:
                sequences = list(SeqIO.parse(compressed_handle, "genbank"))
            assemblies_sequences[assembly['accession']] = sequences
        except (URLError, zlib.error, ValueError) as e:
            console.print(
                f"[bold red]Error for accession {assembly['accession']}: {str(e)}")
            continue

        progress.update(download_task, advance=1, refresh=True)

    return assemblies_sequences


def save_sequences_to_file(assemblies_sequences, progress):
    """Save the downloaded sequences to individual files."""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    assemblies_file_name = f"assemblies_{timestamp}.tsv"
    save_task = progress.add_task(
        "[3/3] Saving sequences to files...", total=len(assemblies_sequences))

    assemblies_info = []

    with open(assemblies_file_name, "w") as assemblies_file:
        assemblies_file.write("AssemblyAccession\tOrganism\tFilePath\n")
        for accession, sequences in assemblies_sequences.items():
            file_path = os.path.abspath(f"sequences/{accession}.gb")
            try:
                os.makedirs("sequences", exist_ok=True)
                with open(file_path, "w") as output_handle:
                    SeqIO.write(sequences, output_handle, "genbank")
                organism = next(
                    iter({seq.annotations["organism"] for seq in sequences}), "Unknown")
                assemblies_file.write(
                    f"{accession}\t{organism}\t{file_path}\n")
                assemblies_info.append((accession, organism, file_path))
            except IOError as e:
                console.print(
                    f"[bold red]Error writing file {file_path}: {str(e)}")
                continue

            progress.update(save_task, advance=1, refresh=True)

    return assemblies_file_name, assemblies_info


def append_to_download_log(assemblies_info):
    """Append the downloaded sequences information to the running log file."""
    with open(DOWNLOAD_LOG_FILE, "a") as log_file:
        if os.stat(DOWNLOAD_LOG_FILE).st_size == 0:
            log_file.write("Timestamp\tAccession\tOrganism\tFilePath\n")

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        for accession, organism, file_path in assemblies_info:
            log_file.write(
                f"{timestamp}\t{accession}\t{organism}\t{file_path}\n")


def display_assemblies_table(assemblies_info):
    """Display the contents of the assemblies information in a rich table."""
    table = Table(title="Assemblies Information")
    table.add_column("Accession", justify="left")
    table.add_column("Organism", justify="left")
    table.add_column("File Path", justify="left")

    for accession, organism, file_path in assemblies_info:
        table.add_row(accession, organism, file_path)

    console.print(table)


def main():
    parser = argparse.ArgumentParser(
        description="Download sequences from NCBI by accession numbers.")
    parser.add_argument("accession_numbers", metavar="ACCESSION",
                        nargs="+", help="An accession number of a sequence to download.")
    parser.add_argument("-s", "--source", choices=["GenBank", "RefSeq"], default="GenBank",
                        help="The source database of the sequences (default: %(default)s).")
    args = parser.parse_args()

    Entrez.email = get_email()

    try:
        with Progress(TextColumn(f"[bold cyan]{{task.description:<{DESCRIPTION_WIDTH}}}", justify="left"), BarColumn(bar_width=None), TextColumn("{task.percentage:>3.0f}%", justify="right"), transient=False) as progress:
            assembly_data = fetch_assembly_data(
                args.accession_numbers, args.source, progress)
            assemblies_sequences = download_assemblies_sequences(
                assembly_data, args.source, progress)
            _, assemblies_info = save_sequences_to_file(
                assemblies_sequences, progress)

        append_to_download_log(assemblies_info)
        console.print(
            f"[bold green]Download complete! {len(assemblies_sequences)} new sequences added to '{DOWNLOAD_LOG_FILE}'.", style="bold green")
        display_assemblies_table(assemblies_info)
    except KeyboardInterrupt:
        console.print("\n[bold red]Interrupted by user. Exiting.")
        sys.exit(1)


if __name__ == "__main__":
    main()

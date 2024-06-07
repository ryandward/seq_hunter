[![DOI](https://zenodo.org/badge/622284759.svg)](https://zenodo.org/badge/latestdoi/622284759)
# Seq Hunter and Seq Dissect NCBI Assembly Sequence Downloader
This repository provides Python scripts to download entire assemblies from the National Center for Biotechnology Information (NCBI) database using accession numbers and to dissect downloaded GenBank files into various formats. The following guide will help you set up and use the scripts effectively.

## Prerequisites
Before you start, ensure you have the following installed:
- Python 3.6 or higher

## Installation and Setup
Clone this repository to your local machine:
```bash
git clone https://github.com/ryandward/seq_hunter.git
cd seq_hunter
```

Create a virtual environment using pipenv:
```bash
pip install pipenv
pipenv install
```

Activate the virtual environment:
```bash
pipenv shell
```

Install the required dependencies:
```bash
pipenv install
```

## Usage

### Seq Hunter
To download nucleotide sequences from NCBI, use the following command:
```bash
python seq_hunter.py <accession_numbers> [--source RefSeq|GenBank]
```
Replace `<accession_numbers>` with one or more space-separated accession numbers (e.g., "NC_000852" "NC_007346"). The `--source` flag is optional and specifies the source of the record (either "RefSeq" or "GenBank"). By default, the script downloads GenBank records.

Example:
```bash
python seq_hunter.py GCF_000005845.2 --source RefSeq
```
The downloaded sequences will be saved as individual GenBank files with the format `<assembly_accession>.gb`. Additionally, a file named `assemblies_<datestamp>.tsv` will be generated, containing a list of accession numbers and corresponding organism names.

### Seq Dissect
To dissect the downloaded GenBank files into various formats (FASTA, BED, FNA, FAA), use the following command:
```bash
python seq_dissect.py <genbank_files>
```
Replace `<genbank_files>` with one or more space-separated GenBank file paths (e.g., "sequences/NC_000852.gb").

Example:
```bash
python seq_dissect.py sequences/GCF_000005845.2.gb
```
This will generate the following files in the `sequences/` directory:
- `<base_name>.fasta`: FASTA format file
- `<base_name>.bed`: BED format file
- `<base_name>.fna`: FNA format file
- `<base_name>.faa`: FAA format file

## Troubleshooting
If you encounter issues, ensure you have the latest version of Python and the required dependencies installed. If problems persist, raise an issue on the GitHub repository.

## License
This project is licensed under the MIT License. Please see the LICENSE file for more details.

seq_hunter NCBI Nucleotide Sequence Downloader
===================================

This repository provides a Python script that allows you to download nucleotide sequences from the National Center for Biotechnology Information (NCBI) database by providing accession numbers. The script is user-friendly, and the following guide will help you set up and use the script effectively.

Prerequisites
-------------

Before you start, please make sure you have the following installed:

- Python 3.6 or higher: Download the latest version of Python from the official website.

Installation and Setup
----------------------

Clone this repository to your local machine:

```bash
git clone https://github.com/yourusername/seq_hunter.git
```

Change to the cloned directory:

```bash
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

Usage
-----

To download nucleotide sequences from NCBI, run the following command:

```bash
python seq_hunter.py <accession_numbers> [--source RefSeq|GenBank]
```

Replace `<accession_numbers>` with one or more space-separated accession numbers (e.g., "NC_000852" "NC_007346"). The `--source` flag is optional and can be used to specify the source of the record (either "RefSeq" or "GenBank"). By default, the script will download GenBank records.

Example:

```bash
python seq_hunter.py NC_000852 NC_007346 --source GenBank
```

Once the script finishes executing, you'll find the downloaded sequences in the current directory, saved as individual GenBank files with the format `<assembly_accession>.gb`. Additionally, a file named `assemblies_<datestamp>.tsv` will be generated, containing a list of the accession numbers and the corresponding organism names.

Troubleshooting
---------------

If you encounter any issues, please make sure you have the latest version of Python and the required dependencies installed. If the problem persists, please raise an issue on the GitHub repository, and we'll help you resolve it.

License
-------

This project is licensed under the MIT License. Please see the LICENSE file for more details.

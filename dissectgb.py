#!/usr/bin/env python3
import sys
import argparse
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
from pathlib import Path

def read_genbank_file(file_path):
    if file_path.suffix == '.gz':
        with gzip.open(file_path, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'genbank'))
    else:
        with file_path.open('r') as handle:
            records = list(SeqIO.parse(handle, 'genbank'))
    return records

# this function extracts all genomic sequences (i.e., chromosomes) 
def write_fasta(records, output_fasta):
    fasta_records = []
    for record in records:
        fasta_records.append(SeqRecord(record.seq, id=record.id, description=record.description))
    SeqIO.write(fasta_records, output_fasta, 'fasta')

# this function finds all bed with a locus tag that are a CDS
def write_bed(records, output_bed):
    with open(output_bed, 'w') as bed_file:
        for record in records:
            for feature in record.features:
                locus_tag = feature.qualifiers.get('locus_tag', [''])[0].strip()
                gene_name = feature.qualifiers.get('gene', [''])[0].strip()
                strand = '+' if feature.location.strand > 0 else '-'
                feature_type = feature.type.strip()
                product = feature.qualifiers.get('product', [''])[0].strip()

                bed_line = [record.id.strip(), str(feature.location.start), str(feature.location.end),
                            locus_tag, gene_name, strand, feature_type, product]
                bed_file.write('\t'.join(bed_line) + '\n')

# this function finds all fna with a locus tag that are a CDS
def write_fna(records, output_fna):
    with open(output_fna, 'w') as fna_file:
        for record in records:
            for feature in record.features:
                locus_tag = feature.qualifiers.get('locus_tag', [''])[0].strip()
                gene_name = feature.qualifiers.get('gene', [''])[0].strip()
                strand = '+' if feature.location.strand > 0 else '-'
                feature_type = feature.type.strip()
                product = feature.qualifiers.get('product', [''])[0].strip()

                if feature_type == 'CDS':
                    fna_line = [record.id.strip(), str(feature.location.start), str(feature.location.end),
                                locus_tag, gene_name, strand, feature_type, product]
                    fna_file.write('>' + locus_tag + '\n')
                    fna_file.write(str(feature.extract(record.seq)) + '\n')

# this function finds all faa with a locus tag that are a CDS
def write_faa(records, output_faa):
    with open(output_faa, 'w') as faa_file:
        for record in records:
            for feature in record.features:
                locus_tag = feature.qualifiers.get('locus_tag', [''])[0].strip()
                gene_name = feature.qualifiers.get('gene', [''])[0].strip()
                strand = '+' if feature.location.strand > 0 else '-'
                feature_type = feature.type.strip()
                product = feature.qualifiers.get('product', [''])[0].strip()

                if feature_type == 'CDS':
                    faa_line = [record.id.strip(), str(feature.location.start), str(feature.location.end),
                                locus_tag, gene_name, strand, feature_type, product]
                    faa_file.write('>' + locus_tag + '\n')
                    faa_file.write(str(feature.qualifiers.get('translation', [''])[0].strip()) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Convert GenBank files to FASTA and BED formats')
    parser.add_argument('input_files', metavar='input_files', nargs='+', help='Input GenBank (gb or gb.gz) files')
    parser.add_argument('--name', type=str, required=True, help='Output file prefix')
    args = parser.parse_args()

    all_records = []
    for input_file in args.input_files:
        all_records.extend(read_genbank_file(Path(input_file)))

    output_fasta = f"{args.name}.fasta"
    output_bed = f"{args.name}.bed"
    output_fna = f"{args.name}.fna"
    output_faa = f"{args.name}.faa"

    write_fasta(all_records, output_fasta)
    write_bed(all_records, output_bed)
    write_fna(all_records, output_fna)
    write_faa(all_records, output_faa)

if __name__ == '__main__':
    main()

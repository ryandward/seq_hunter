import sys
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

def write_fasta(records, output_fasta):
    fasta_records = []
    for record in records:
        fasta_records.append(SeqRecord(record.seq, id=record.id, description=record.description))
    SeqIO.write(fasta_records, output_fasta, 'fasta')

def write_bed(records, output_bed):
    with open(output_bed, 'w') as bed_file:
        for record in records:
            for feature in record.features:
                locus_tag = feature.qualifiers.get('locus_tag', [''])[0].strip()
                gene_name = feature.qualifiers.get('gene', [''])[0].strip()
                strand = '+' if feature.location.strand > 0 else '-'
                feature_type = feature.type.strip()
                product = feature.qualifiers.get('product', [''])[0].strip()

                if isinstance(feature.location, FeatureLocation):
                    start, end = feature.location.start, feature.location.end
                    bed_line = [record.id.strip(), str(start), str(end), locus_tag, gene_name, strand, feature_type, product]
                    bed_file.write('\t'.join(bed_line) + '\n')
                else:  # CompoundLocation
                    for part in feature.location.parts:
                        start, end = part.start, part.end
                        bed_line = [record.id.strip(), str(start), str(end), locus_tag, gene_name, strand, feature_type, product]
                        bed_file.write('\t'.join(bed_line) + '\n')

def write_fna(records, output_fna):
    with open(output_fna, 'w') as fna_file:
        for record in records:
            for feature in record.features:
                locus_tag = feature.qualifiers.get('locus_tag', [''])[0].strip()
                gene_name = feature.qualifiers.get('gene', [''])[0].strip()
                product = feature.qualifiers.get('product', [''])[0].strip()
                feature_type = feature.type.strip()

                if feature_type == 'CDS':
                    header = f">{locus_tag}|{gene_name}|{gene_name.upper()}{product}"
                    fna_file.write(header + '\n')
                    fna_file.write(str(feature.extract(record.seq)) + '\n')

def write_faa(records, output_faa):
    with open(output_faa, 'w') as faa_file:
        for record in records:
            for feature in record.features:
                locus_tag = feature.qualifiers.get('locus_tag', [''])[0].strip()
                gene_name = feature.qualifiers.get('gene', [''])[0].strip()
                product = feature.qualifiers.get('product', [''])[0].strip()
                feature_type = feature.type.strip()

                if feature_type == 'CDS':
                    header = f">{locus_tag}|{gene_name}|{product}"
                    faa_file.write(header + '\n')
                    faa_file.write(str(feature.qualifiers.get('translation', [''])[0].strip()) + '\n')

def main():
    for input_file in sys.argv[1:]:
        all_records = read_genbank_file(Path(input_file))

        # Set output file names based on the input file name
        input_file_name = Path(input_file).stem
        output_fasta = f"sequences/{input_file_name}.fasta"
        output_bed = f"sequences/{input_file_name}.bed"
        output_fna = f"sequences/{input_file_name}.fna"
        output_faa = f"sequences/{input_file_name}.faa"

        # Write the output files for the current input file
        write_fasta(all_records, output_fasta)
        write_bed(all_records, output_bed)
        write_fna(all_records, output_fna)
        write_faa(all_records, output_faa)

if __name__ == '__main__':
    main()

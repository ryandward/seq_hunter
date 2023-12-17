import sys
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
from pathlib import Path


def read_genbank_file(file_path):
    mode = 'rt' if file_path.suffix == '.gz' else 'r'
    with gzip.open(file_path, mode) if file_path.suffix == '.gz' else file_path.open(mode) as handle:
        return list(SeqIO.parse(handle, 'genbank'))


def write_fasta(records, output_fasta):
    fasta_records = [SeqRecord(
        record.seq, id=record.id, description=record.description) for record in records]
    SeqIO.write(fasta_records, output_fasta, 'fasta')


def format_feature(feature, record_id):
    locus_tag = feature.qualifiers.get('locus_tag', [''])[0].strip()
    gene_name = feature.qualifiers.get('gene', [''])[0].strip()
    strand = '+' if feature.location.strand > 0 else '-'
    feature_type = feature.type.strip()
    product = feature.qualifiers.get('product', [''])[0].strip()
    start, end = feature.location.start, feature.location.end
    return [record_id.strip(), str(start), str(end), locus_tag, gene_name, strand, feature_type, product]


def write_bed(records, output_bed):
    with open(output_bed, 'w') as bed_file:
        for record in records:
            sorted_features = sorted(
                record.features, key=lambda f: f.location.start)
            for feature in sorted_features:
                if feature.type == "CDS" and not feature.qualifiers.get('translation'):
                    continue  # Skip CDS features with empty protein sequences
                if isinstance(feature.location, FeatureLocation):
                    bed_line = format_feature(feature, record.id)
                    bed_file.write('\t'.join(bed_line) + '\n')
                else:  # CompoundLocation
                    for part in feature.location.parts:
                        bed_line = format_feature(feature, record.id)
                        bed_line[1] = str(part.start)  # Start position
                        bed_line[2] = str(part.end)    # End position
                        bed_file.write('\t'.join(bed_line) + '\n')


def format_header(locus_tag, gene_name, product):
    gene_name = gene_name + ' ' if gene_name else ''
    return f">{locus_tag}|{gene_name}{product}"

def write_feature_sequence(output_file, records, feature_type, sequence_getter):
    with open(output_file, 'w') as file:
        for record in records:
            for feature in record.features:
                if feature.type.strip() == feature_type and 'translation' in feature.qualifiers and feature.qualifiers.get('translation'):
                    locus_tag = feature.qualifiers.get('locus_tag', [''])[0].strip()
                    gene_name = feature.qualifiers.get('gene', [''])[0].strip()
                    product = feature.qualifiers.get('product', [''])[0].strip()
                    header = format_header(locus_tag, gene_name, product)
                    file.write(header + '\n')
                    file.write(sequence_getter(feature, record) + '\n')

def write_fna(records, output_fna):
    write_feature_sequence(output_fna, records, 'CDS', lambda feature, record: str(feature.extract(record.seq)))

def write_faa(records, output_faa):
    write_feature_sequence(output_faa, records, 'CDS', lambda feature, _: feature.qualifiers.get('translation', [''])[0].strip())


def main():
    for input_file in sys.argv[1:]:
        all_records = read_genbank_file(Path(input_file))
        base_name = Path(input_file).stem

        write_fasta(all_records, f"sequences/{base_name}.fasta")
        write_bed(all_records, f"sequences/{base_name}.bed")
        write_fna(all_records, f"sequences/{base_name}.fna")
        write_faa(all_records, f"sequences/{base_name}.faa")


if __name__ == '__main__':
    main()

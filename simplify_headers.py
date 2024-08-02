import argparse
from Bio import SeqIO

def simplify_headers(input_file, output_file, mapping_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(mapping_file, 'w') as mapfile:
        for i, record in enumerate(SeqIO.parse(infile, 'fasta'), start=1):
            original_id = record.id
            new_id = "seq{}".format(i)
            record.id = new_id
            record.description = ""
            SeqIO.write(record, outfile, 'fasta')
            mapfile.write("{}\t{}\n".format(new_id, original_id))

def main():
    parser = argparse.ArgumentParser(description="Simplify FASTA headers and create a mapping file.")
    parser.add_argument('input_fasta', help="Input FASTA file")
    parser.add_argument('output_fasta', help="Output FASTA file with simplified headers")
    parser.add_argument('mapping_file', help="Output mapping file")

    args = parser.parse_args()

    simplify_headers(args.input_fasta, args.output_fasta, args.mapping_file)

if __name__ == "__main__":
    main()


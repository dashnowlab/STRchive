import csv
import json
import sys
import argparse
sys.path.append('../')
import reftogeneorientation

def main(in_csv, out_csv, out_json):
    # This will create the gene orientation columns from ref orientation for each motif type
    # Reference/benign, pathogenic, and unknown
    reftogeneorientation.process_csv(in_csv, out_csv)

    # Read CSV file, created by process_csv
    with open(out_csv, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        data = list(csv_reader)

    # Write JSON file
    with open(out_json, 'w') as json_file:
        json_file.write(json.dumps(data, indent=4))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process CSV and convert to JSON')
    parser.add_argument('--in_csv', default='data/STR-disease-loci.csv', help='Input CSV file path')
    parser.add_argument('--out_csv', default='data/STR-disease-loci.processed.csv', help='Output CSV file path')
    parser.add_argument('--out_json', default='data/STR-disease-loci.processed.json', help='Output JSON file path')

    args = parser.parse_args()
    main(args.in_csv, args.out_csv, args.out_json)

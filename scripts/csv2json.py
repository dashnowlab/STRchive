import csv
import json
import sys
sys.path.append('../')
import reftogeneorientation

def main():
    #Here we are defining the original CSV as well as the output csv/json
    in_csv = 'data/STR-disease-loci.csv'
    out_json = 'data/STR-disease-loci.processed.json'
    out_csv = 'data/STR-disease-loci.processed.csv'

    #This will create the gene orientation columns from ref orientation for each motif type
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
    main()
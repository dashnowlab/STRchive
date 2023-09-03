import csv
import json

def main():

    in_csv = 'data/STR-disease-loci.csv'
    out_json = 'data/STR-disease-loci.json'

    # Read CSV file
    with open(in_csv, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        data = list(csv_reader)

    # Write JSON file
    with open(out_json, 'w') as json_file:
        json_file.write(json.dumps(data, indent=4))

if __name__ == '__main__':
    main()
# Overwrite the CSV or JSON file, whichever is older

import csv
import json
import sys
import argparse
import doctest
import os
import pandas as pd
# in case person runs from data instead of scripts
sys.path.append('../')
import reftogeneorientation

def main(csv_fname, json_fname):

    # Check if the CSV file is newer than the JSON file
    if os.path.getmtime(csv_fname) > os.path.getmtime(json_fname):
        sys.stderr.write('CSV is newer, reading CSV and writing JSON\n')
        df = pd.read_csv(csv_fname)
        # Create/update gene orientation columns from ref orientation for each motif type:
        # Reference/benign, pathogenic, and unknown
        df = reftogeneorientation.process_df(df)
        # Sort by gene name
        df.sort_values(by='gene', inplace=True)
        # Write JSON file
        with open(json_fname, 'w') as json_file:
            json_string = df.to_json(orient='records', indent=4)
            json_file.write(json_string)

    else:
        sys.stderr.write('JSON is newer, reading JSON and writing CSV\n')
        # Check if file exists
        if not os.path.exists(json_fname):
            sys.stderr.write(f'Error: {json_fname} does not exist\n')
            sys.exit(1)
        df = pd.read_json(json_fname, orient='records')
        # Create/update gene orientation columns from ref orientation for each motif type:
        # Reference/benign, pathogenic, and unknown
        df = reftogeneorientation.process_df(df)
        # Sort by gene name
        df.sort_values(by='gene', inplace=True)
        # Write CSV file
        with open(csv_fname, 'w') as csv_file:
            df.to_csv(csv_file, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read CSV and JSON versions of the database and overwrite the older one to match the newer one.')
    parser.add_argument('--csv', default='data/STRchive-database.csv', help='CSV file path')
    parser.add_argument('--json', default='data/STRchive-database.json', help='JSON file path')
    doctest.testmod()
    args = parser.parse_args()
    main(args.csv, args.json)

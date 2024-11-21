# Overwrite the CSV or JSON file, whichever is older

import sys
import argparse
import doctest
import os
import pandas as pd
import time

# in case person runs from data instead of scripts
sys.path.append('../')
import reftogeneorientation

def main(csv_fname, json_fname, pause):

    # Check if the CSV file is newer than the JSON file
    if os.path.getmtime(csv_fname) > os.path.getmtime(json_fname):
        sys.stderr.write('CSV is newer, reading CSV and writing JSON\n')
        sys.stderr.write(f'WARNING: overwriting {json_fname} in {pause} seconds\n')
        sys.stderr.write('Press Ctrl+C to cancel\n')
        time.sleep(pause)

        df = pd.read_csv(csv_fname)
        # Create/update gene orientation columns from ref orientation for each motif type:
        # Reference/benign, pathogenic, and unknown
        df = reftogeneorientation.process_df(df)
        # Sort by gene name
        df.sort_values(by=['gene', 'id'], inplace=True)
        # Write JSON file
        with open(json_fname, 'w') as json_file:
            json_string = df.to_json(orient='records', indent=4)
            json_file.write(json_string)

    else:
        sys.stderr.write('JSON is newer, reading JSON and writing CSV\n')
        sys.stderr.write(f'WARNING: overwriting {csv_fname} in {pause} seconds\n')
        sys.stderr.write('Press Ctrl+C to cancel\n')
        time.sleep(pause)

        # Check if file exists
        if not os.path.exists(json_fname):
            sys.stderr.write(f'Error: {json_fname} does not exist\n')
            sys.exit(1)
        df = pd.read_json(json_fname, orient='records')
        # Create/update gene orientation columns from ref orientation for each motif type:
        # Reference/benign, pathogenic, and unknown
        df = reftogeneorientation.process_df(df)
        # Sort by gene name
        df.sort_values(by=['gene', 'id'], inplace=True)
        # Write CSV file
        with open(csv_fname, 'w') as csv_file:
            df.to_csv(csv_file, index=False)
    
    # Update the TRGT repeat definitions by running make_catalog.py
    # ref_genomes = ['hg19', 'hg38', 'T2T']
    # for genome in ref_genomes:
    #     os.system(f'python scripts/make_catalog.py -g {genome} -f TRGT data/STRchive-database.csv data/{genome}.STRchive-disease-loci.TRGT.bed')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read CSV and JSON versions of the database and overwrite the older one to match the newer one.')
    parser.add_argument('--csv', default='data/STRchive-database.csv', help='CSV file path')
    parser.add_argument('--json', default='data/STRchive-database.json', help='JSON file path')
    parser.add_argument('--pause', default=5, help='Pause time in seconds before overwriting the file')
    doctest.testmod()
    args = parser.parse_args()
    main(args.csv, args.json, args.pause)

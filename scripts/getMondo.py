# Download the latest Mondo disease descriptions using their rest API
# Check for updates

import requests # Need to add to environment.yml
import sys
import json

mondo_url_prefix = "https://www.ebi.ac.uk/ols4/api/ontologies/mondo/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FMONDO_"
strchive_json = "data/STRchive-database.json"

def getMondoIDs(strchive_json):
    with open(strchive_json, 'r') as f:
        data = json.load(f)
    MondoIDs = []
    for entry in data:
        MondoIDs.append(entry['Mondo'])
    return MondoIDs

def getMondoDesc(MondoID):
    try:
        response = requests.get(mondo_url_prefix + MondoID)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        sys.err('Error: ' + str(e) + '\n')
        return None
    description = getDiseaseDescriptions(response.json())
    return description

def getDiseaseDescriptions(mondo_json):
    description_list = mondo_json['description']
    if len(description_list) == 0:
        return None
    return description_list[0]

def main():
    # Iterate through the STRchive database and get the Mondo IDs
    with open(strchive_json, 'r') as f:
        data = json.load(f)
    for entry in data:
        MondoIDs = entry['Mondo']
        descriptions = []
        if MondoIDs is not None:
            for MondoID in MondoIDs.split(';'):
                description = getMondoDesc(MondoID.strip())
                if description is not None:
                    descriptions.append(description)
            
        # Add/update the descriptions in the STRchive database
        #print(entry['gene'], descriptions)
        if len(descriptions) > 0:
            entry['disease_description'] = '; '.join(descriptions)
        else:
            entry['disease_description'] = None


    # Write the updated STRchive database
    with open(strchive_json, 'w') as f:
        # Write out json with indenting
        json.dump(data, f, indent=4, separators=(',', ':'), ensure_ascii=True)


if __name__ == '__main__':
    main()
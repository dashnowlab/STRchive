#!/bin/bash

# Input and output files
input_file=$1
output_file=$2

# Start JSON array
echo "[" > "$output_file"
first=true

# Iterate over each citation
while IFS= read -r citation; do
    if [ "$first" = true ]; then
        first=false
    else
        echo "," >> "$output_file"  # Separate JSON objects with a comma
    fi

    # Try to fetch the JSON from manubot
    citation_json=$(manubot cite "$citation" 2>/dev/null)
    if [ $? -eq 0 ]; then
        # If successful, add the fetched JSON
        echo "$citation_json" >> "$output_file"
    else
        # If unsuccessful, add an entry with the original citation and an error message
        echo "{ \"citation\": \"$citation\", \"error\": \"Unable to parse\" }" >> "$output_file"
    fi
done < <(jq -r '.[]' "$input_file")

# End JSON array
echo "]" >> "$output_file"

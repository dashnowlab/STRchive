#!/bin/bash

# Iterate over each citation passed as an argument
echo "["  # Start the JSON array
first=true
for citation in "$@"; do
    if [ "$first" = true ]; then
        first=false
    else
        echo ","  # Separate JSON objects with a comma
    fi
    # Fetch the JSON from manubot
    manubot cite "$citation"
done
echo "]"  # End the JSON array

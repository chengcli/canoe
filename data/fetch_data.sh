#! /bin/bash

# Download hitran data from Dropbox
FILE="HITRAN2012.par"

# Dropbox link
DROPBOX_LINK="https://www.dropbox.com/s/5tvhdll2yrpnxue/HITRAN2012.par.tar.gz?dl=0"

# Read expected SHA256 from the file
EXPECTED_SHA256=$(grep "$FILE" checksums.txt | awk '{print $1}')

# Check if the SHA256 has been found
if [ -z "$EXPECTED_SHA256" ]; then
    echo "Cannot find SHA256 for $FILE in checksums.txt"
    exit 1
fi

# Check if the file exists
if [ -f "$FILE" ]; then
    echo "$FILE exists. Checking SHA256..."
    
    # Check the operating system and calculate SHA256 of the existing file
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        SHA256=$(sha256sum $FILE | cut -d ' ' -f1)
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        SHA256=$(shasum -a 256 $FILE | cut -d ' ' -f1)
    else
        echo "Unsupported operating system"
        exit 1
    fi

    # Compare the calculated SHA256 with the expected one
    if [ "$SHA256" == "$EXPECTED_SHA256" ]; then
        echo "SHA256 is correct. No need to download the file."
    else
        echo "SHA256 is not correct. Need to download the file..."
    fi
else
    echo "$FILE does not exist. Need to download the file..."
fi

# Download and extract the file if it does not exist or its SHA256 is incorrect
if [ ! -f "$FILE" ] || [ "$SHA256" != "$EXPECTED_SHA256" ]; then
    wget -q --show-progress -O HITRAN2012.par.tar.gz $DROPBOX_LINK
    if tar -xzvf HITRAN2012.par.tar.gz; then
        echo "Successfully extracted HITRAN2012.par.tar.gz"
        rm HITRAN2012.par.tar.gz
        echo "Removed HITRAN2012.par.tar.gz"
    else
        echo "Failed to extract HITRAN2012.par.tar.gz"
        exit 1
    fi
fi

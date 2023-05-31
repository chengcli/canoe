#! /bin/bash

# Download cia data from Dropbox
TAR_FILE="H2-He-cia"
DATA_DIR=$1/

# Dropbox link
DROPBOX_LINK="https://www.dropbox.com/s/czbzu9hglxty4rj/H2-He-cia.tar.gz?dl=0"

# Get all files from the checksum file
FILES=(H2-H2-eq.orton.txt H2-H2-eq.xiz.txt \
    H2-H2-nm.orton.txt H2-H2-nm.xiz.txt \
    H2-He-eq.orton.txt H2-He-eq.xiz.txt \
    H2-He-nm.orton.txt H2-He-nm.xiz.txt)

# Function to calculate SHA256 based on the OS
function calculate_sha256 {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo $(sha256sum $1 | cut -d ' ' -f1)
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        echo $(shasum -a 256 $1 | cut -d ' ' -f1)
    else
        echo "Unsupported operating system"
        exit 1
    fi
}

# Check all files
DOWNLOAD=FALSE
for FILE in ${FILES[@]}; do
    # Read expected SHA256 from the file
    EXPECTED_SHA256=$(grep "$FILE" ${DATA_DIR}checksums.txt | awk '{print $1}')

    # Check if the SHA256 has been found
    if [ -z "$EXPECTED_SHA256" ]; then
        echo "Cannot find SHA256 for $FILE in ${DATA_DIR}checksums.txt"
        exit 1
    fi

    # Check if the file exists
    if [ -f "${DATA_DIR}${FILE}" ]; then
        echo "$FILE exists. Checking SHA256..."
        
        SHA256=$(calculate_sha256 ${DATA_DIR}${FILE})
    
        # Compare the calculated SHA256 with the expected one
        if [ "$SHA256" != "$EXPECTED_SHA256" ]; then
            echo "SHA256 for $FILE is not correct. Need to download and extract the tar file again..."
            DOWNLOAD=TRUE
            break
        else
            echo "SHA256 for $FILE is correct. No need to download the file"
        fi
    else
        echo "$FILE does not exist. Need to download and extract the tar file..."
        DOWNLOAD=TRUE
        break
    fi
done

# Check if any file is missing or has incorrect SHA256
if [ "$DOWNLOAD" = "TRUE" ]; then
    wget -q --show-progress -O ${DATA_DIR}${TAR_FILE}.tar.gz $DROPBOX_LINK
    if tar -xzvf ${DATA_DIR}${TAR_FILE}.tar.gz -C ${DATA_DIR}; then
        echo "Successfully extracted ${TAR_FILE}.tar.gz"
        rm ${DATA_DIR}${TAR_FILE}.tar.gz
        echo "Removed ${TAR_FILE}.tar.gz"
    else
        echo "Failed to extract ${TAR_FILE}.tar.gz"
        exit 1
    fi
fi

#!/bin/bash

# Set the base directory
BASE_DIR=$(dirname "$(dirname "$(realpath $0)")")

# Directories and files to remove
DIRS_TO_REMOVE=(
    "${BASE_DIR}/pipelines/cache"
    "${BASE_DIR}/pipelines/temp"
    "${BASE_DIR}/pipelines/temp_adata"
)

# Remove directories and files
echo "Removing temporary directories and files:"

for DIR in "${DIRS_TO_REMOVE[@]}"; do
    if [ -d "$DIR" ]; then
        rm -rf "$DIR"
        echo "Removed: $DIR"
    else
        echo "Directory does not exist: $DIR"
    fi
done

echo ""

# Directories in results to remove
RES_DIRS_TO_REMOVE=$(find "${BASE_DIR}/results" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;)

# Remove directories in results dir
echo "Removing directories in results directory:"

for DIR in $RES_DIRS_TO_REMOVE; do
    DIR_PATH="${BASE_DIR}/results/$DIR"

    if [ -d "$DIR_PATH" ]; then
        rm -rf "$DIR_PATH"
        echo "Removed directory: $DIR_PATH"
    fi
done

echo ""
echo "Cleanup complete."
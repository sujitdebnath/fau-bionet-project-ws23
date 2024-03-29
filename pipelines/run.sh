#!/bin/bash

# Find the project directory and set the disease and datset id
PROJECT_DIR=$(dirname "$(dirname "$(realpath $0)")")
disease_id="mpn"
dataset_id="dataset1"

# Run the data handler script with required arguments
echo "----------------------------------- Data Extraction Started -----------------------------------"
python3 "${PROJECT_DIR}/pipelines/services/data_handler.py" \
    --disease_id "${disease_id}" \
    --dataset_id "${dataset_id}"
echo "----------------------------------- Data Extraction Ended -----------------------------------"

echo ""
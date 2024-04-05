#!/bin/bash

# Find the project directory
PROJECT_DIR=$(dirname "$(dirname "$(realpath $0)")")

# Find all the diseases in the dataset dir
diseases=$(find "${PROJECT_DIR}/dataset" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;)

for disease_id in $diseases; do
    # Find all the dataset available for a specific disease
    datasets=$(find "${PROJECT_DIR}/dataset/${disease_id}" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;)

    # run whole pipeline for all dataset of a specific disease
    for dataset_id in $datasets; do
        echo "---------------------------- Pipelines started for disease: $disease_id, data: $dataset_id ----------------------------"

        # Run the data handler script with required arguments
        echo "----------------------------------- Data Handler Started -----------------------------------"
        python3 "${PROJECT_DIR}/pipelines/services/adata_handler.py" \
            --disease_id "${disease_id}" \
            --dataset_id "${dataset_id}"
        echo "----------------------------------- Data Handler Ended -----------------------------------"

        echo ""

        # Run the data preprocessing script with required arguments
        echo "----------------------------------- Data Preprocessing Started -----------------------------------"
        python3 "${PROJECT_DIR}/pipelines/services/adata_preprocessor.py" \
            --disease_id "${disease_id}" \
            --dataset_id "${dataset_id}" \
            --mito_perc 0.05 \
            --n_umis 500 \
            --detected_genes 250 \
            --preprocess_mode 'shiftlog|pearson' \
            --n_hvgs 2000 \
            --n_pcs 50 \
            --n_neighbors 15 \
            --cluster_mode 'leiden'
        echo "----------------------------------- Data Preprocessing Ended -----------------------------------"

        echo ""

        # Run the cell type annotation script with required arguments
        echo "----------------------------------- Cell-type Annotation Started -----------------------------------"
        python3 "${PROJECT_DIR}/pipelines/services/cell_type_annotation.py" \
            --disease_id "${disease_id}" \
            --dataset_id "${dataset_id}" \
            --anno_methods "scsa" "metatime" \
            --cellmarker \
            --panglaodb \
            --foldchange 1.5 \
            --pvalue 0.01 \
            --celltype 'normal' \
            --clustertype 'leiden' \
            --rank_rep
        echo "----------------------------------- Cell-type Annotation Ended -----------------------------------"

        echo ""

        # Run the differential gene expression analysis script with required arguments
        echo "----------------------------------- DGE Analysis Started -----------------------------------"
        python3 "${PROJECT_DIR}/pipelines/services/diff_gene_exp_analysis.py" \
            --disease_id "${disease_id}" \
            --dataset_id "${dataset_id}" \
            --dge_methods "t-test" "wilcoxon" "logreg" "t-test_overestim_var"
        echo "----------------------------------- DGE Analysis Ended -----------------------------------"

        echo "---------------------------- Pipelines ended for disease: $disease_id, data: $dataset_id ----------------------------\n"
    done
done

# Unnecessary directories and files to remove
DIRS_TO_REMOVE=(
    "${PROJECT_DIR}/pipelines/cache"
    "${PROJECT_DIR}/pipelines/temp"
    "${PROJECT_DIR}/pipelines/temp_adata"
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

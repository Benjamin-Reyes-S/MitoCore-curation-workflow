#!/usr/bin/env bash
# Downloads the input models and external database files the curation pipeline
# needs, skipping anything that is already present (e.g. on a mounted volume).
set -e

INPUT_MODELS_DIR="Input_Models"
FILES_DATABASES_DIR="Files_Databases"

mkdir -p "$INPUT_MODELS_DIR" "$FILES_DATABASES_DIR"

# fetch <url> <destination>
# Skips the download if the destination already exists and is non-empty.
# Downloads to a temp file first so a failed/interrupted transfer never leaves
# a partial file at the destination.
fetch() {
    local url="$1"
    local dest="$2"

    if [ -s "$dest" ]; then
        echo "skip (already present): $dest"
        return 0
    fi

    echo "downloading: $dest"
    echo "  from: $url"
    curl -fSL --retry 3 --retry-delay 5 -o "${dest}.tmp" "$url"
    mv "${dest}.tmp" "$dest"
}

# ---- models (mdoa-group/curating-mitocore, the originals) ----
fetch "https://gitlab.com/mdoa-group/curating-mitocore/-/raw/main/model_curation/Input_Models/MitoCore_Original_2017.xml" \
    "$INPUT_MODELS_DIR/MitoCore_Original_2017.xml"

fetch "https://gitlab.com/mdoa-group/curating-mitocore/-/raw/main/model_curation/Input_Models/MitoMAMMAL_08_25.xml" \
    "$INPUT_MODELS_DIR/MitoMAMMAL_08.25.xml"

fetch "https://gitlab.com/mdoa-group/curating-mitocore/-/raw/main/model_curation/Input_Models/Human-GEM.xml" \
    "$INPUT_MODELS_DIR/Human-GEM.xml"

# ---- BiGG metabolite bulk mapping ----
# Note: this is the *metabolites* namespace dump (bigg_models_metabolites.txt),
# with the 'universal_bigg_id'/'database_links' columns PreliminaryCuration.py
# reads. bigg_models_reactions.txt is a different file and won't work here.
fetch "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt" \
    "$FILES_DATABASES_DIR/bigg_models_metabolites.txt"

# ---- MetaNetX (latest MNXref release) ----
fetch "https://www.metanetx.org/ftp/latest/chem_xref.tsv" \
    "$FILES_DATABASES_DIR/MNX_chem_xref.tsv"

fetch "https://www.metanetx.org/ftp/latest/reac_xref.tsv" \
    "$FILES_DATABASES_DIR/MNX_reac_xref.tsv"

fetch "https://www.metanetx.org/ftp/latest/reac_prop.tsv" \
    "$FILES_DATABASES_DIR/reac_prop.tsv"

# ---- Human-GEM id-mapping tables ----
fetch "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/reactions.tsv" \
    "$FILES_DATABASES_DIR/Human1_reactions.tsv"

fetch "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/metabolites.tsv" \
    "$FILES_DATABASES_DIR/Human1_metabolites.tsv"

fetch "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/genes.tsv" \
    "$FILES_DATABASES_DIR/Human1_genes.tsv"

echo "All input models and database files are present."

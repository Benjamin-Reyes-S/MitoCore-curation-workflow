#!/usr/bin/env bash
# Container entrypoint: fetch the mounted-volume data (if not already there),
# then run the curation pipeline.
set -e

./download_data.sh
./config.sh

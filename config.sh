
#!/usr/bin/env bash
echo "Starting Curation"

set -e # makes stop on error inmeditely 

echo "Starting Curation"

if [ ! -d ".venv" ]; then #if no virtual env
    python3 -m venv .venv #create virtual env
    . .venv/bin/activate # . means source this file in the current shell
    pip install -r requirements.txt #install python req

else
    . .venv/bin/activate #else activate virtual env
fi

export PYTHONPATH="$(pwd)"
python3 -m OriginalParse.ParsingOriginalModel

python3 -m PreliminaryCuration.PreliminaryCuration

python3 -m MitoMammalToMitoCore.MitoMammalToMitoCore

#python3 -m AlignmentToHuman1.AlignmentToHuman1


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
python3 -m Original_Parse.Parsing_Original_Model


#MitoCore preliminary curation
#docker compose -f folde/docker-compose.yml up -d --wait
#python -m Preliminary_Curation.method

#MitoCore MitoMammal parse 
#docker compose -f folder/docker-compose.yml up -d --wait
#python -m MitoMammal_Parse.method

#MitoCore MitoMammal parse 
#docker compose -f folder/docker-compose.yml up -d --wait
#python -m MitoMammal_Curation.method

#MitoCore alignment to Human1
#docker compose -f folder/docker-compose.yml up -d --wait
#python -m Human1_Alignment.method
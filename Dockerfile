FROM python:3.11-slim

WORKDIR /app

# download_data.sh needs curl; libexpat1 is the shared lib python-libsbml's
# compiled extension links against; git is required by memote's GitPython
# dependency. None of these ship with the slim base image.
RUN apt-get update \
    && apt-get install -y --no-install-recommends curl libexpat1 git \
    && rm -rf /var/lib/apt/lists/*

# Install dependencies first so this layer is cached across code-only changes.
COPY requirements.txt .
RUN python3 -m venv .venv \
    && . .venv/bin/activate \
    && pip install --no-cache-dir -r requirements.txt

# Pipeline code only -- input/output models and database files are supplied
# at runtime via mounted volumes (see docker-compose.yml), never baked in.
COPY OriginalParse/ OriginalParse/
COPY PreliminaryCuration/ PreliminaryCuration/
COPY MitoMammalToMitoCore/ MitoMammalToMitoCore/
COPY AlignmentToHuman1/ AlignmentToHuman1/
COPY src/ src/
COPY config.sh download_data.sh docker-entrypoint.sh ./
ENTRYPOINT ["./docker-entrypoint.sh"]

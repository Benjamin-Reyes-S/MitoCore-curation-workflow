FROM python:3.11-slim

WORKDIR /app

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

RUN chmod +x config.sh download_data.sh docker-entrypoint.sh

ENTRYPOINT ["./docker-entrypoint.sh"]

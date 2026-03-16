FROM python:3.8-slim

# Install system dependencies + GIT
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl build-essential libsuitesparse-dev git \
    && rm -rf /var/lib/apt/lists/*

# Install Poetry
ENV POETRY_HOME="/opt/poetry"
ENV PATH="$POETRY_HOME/bin:$PATH"
RUN curl -sSL https://install.python-poetry.org | python3 - --version 1.3.2

WORKDIR /app

# 1. Clone the repository directly into the container
# 2. Initialize submodules (as seen in your YAML)
RUN git clone https://github.com/GeomScale/dingo.git . \
    && git submodule update --init

# 3. Download the C++ libraries (Boost & lp-solve)
RUN wget -O boost.tar.bz2 https://archives.boost.io/release/1.76.0/source/boost_1_76_0.tar.bz2 \
    && tar xjf boost.tar.bz2 && rm boost.tar.bz2
    
RUN wget https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_source.tar.gz \
    && tar xzvf lp_solve_5.5.2.11_source.tar.gz && rm lp_solve_5.5.2.11_source.tar.gz

# 4. Install Python dependencies
RUN poetry config virtualenvs.create false \
    && poetry install --no-interaction \
    && pip3 install numpy scipy jupyter

CMD ["python3"]
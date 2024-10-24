# Skye Goetz (CalPoly) 10/23/24

# Base image
FROM julia:1.9

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    python3-pip \
    python3-dev \
    python3-venv \
    libtiff5-dev \
    libboost-python-dev \
    cmake \
    git \
    libgomp1 \
    gnupg2 \
    software-properties-common \
    && rm -rf /var/lib/apt/lists/*

# Setup Python virtual environment
RUN python3 -m venv /opt/venv \
    && . /opt/venv/bin/activate \
    && pip install --upgrade pip \
    && pip install numpy pandas rdkit

ENV PATH="/opt/venv/bin:$PATH"

COPY orca.tar.xz /opt/orca.tar.xz

# Extract and install ORCA [[ Assumes Download Named orca.tar.xz Is In The Same Directory ]]
WORKDIR /opt
RUN tar -xJf /opt/orca.tar.xz \
    && rm /opt/orca.tar.xz

ENV PATH="/opt/orca:/opt/orca/bin:$PATH"

RUN julia -e 'using Pkg; Pkg.add(["FilePathsBase", "LaTeXStrings", "DataFrames", "PyCall", "Plots", "YAML", "CSV"]); Pkg.precompile()'

WORKDIR /app

COPY . /app

CMD ["julia", "OrcaFunctionalHub.jl"]
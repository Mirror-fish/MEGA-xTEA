# ============================================================
# MEGA-xTEA Docker Image
# Multi-stage build: compile C++ core → slim runtime image
# NOTE: Ubuntu 20.04 (Python 3.8) required for deep-forest
#       compatibility — deep-forest only has PyPI wheels for
#       Python 3.7/3.8/3.9 (last release Sept 2022).
# ============================================================

# ---------- Stage 1: Build ----------
FROM ubuntu:20.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    g++ \
    make \
    autoconf \
    automake \
    libtool \
    pkg-config \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Build htslib 1.19 from source
WORKDIR /tmp
RUN wget -q https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 && \
    tar xjf htslib-1.19.tar.bz2 && \
    cd htslib-1.19 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    ldconfig

# Build samtools 1.19 from source (apt samtools 1.10 on 20.04 lacks `depth -@`)
RUN wget -q https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar xjf samtools-1.19.tar.bz2 && \
    cd samtools-1.19 && \
    ./configure --prefix=/usr/local --with-htslib=/usr/local && \
    make -j$(nproc) && \
    make install

# Copy MEGA-xTEA source
WORKDIR /opt/mega-xtea
COPY cpp/ cpp/
COPY Makefile .

# Patch Makefile to use system htslib instead of external/htslib
RUN sed -i 's|LHTS = $(CURDIR)/external/htslib|LHTS = /usr/local|' Makefile && \
    sed -i 's|-I $(LHTS)|-I /usr/local/include|g' Makefile && \
    sed -i 's|-L $(LHTS)|-L /usr/local/lib|g' Makefile && \
    sed -i 's|-Wl,-rpath=$(LHTS)||g' Makefile

# Compile C++ modules
RUN make -j$(nproc)

# Verify compiled binaries
RUN ls -la cpp/extract_discordant cpp/extract_unmapped cpp/*.so

# ---------- Stage 2: Runtime ----------
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Runtime dependencies (samtools built from source, not apt)
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    build-essential \
    python3-dev \
    bedtools \
    ncbi-blast+ \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4 \
    && rm -rf /var/lib/apt/lists/*

# Install htslib + samtools 1.19 from builder
COPY --from=builder /usr/local/lib/libhts* /usr/local/lib/
COPY --from=builder /usr/local/include/htslib /usr/local/include/htslib
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/samtools
RUN ldconfig

# Install Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir "pip<24.1" "setuptools<58.0" wheel && \
    pip3 install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

# Copy MEGA-xTEA
WORKDIR /opt/mega-xtea
COPY --from=builder /opt/mega-xtea/cpp/extract_discordant cpp/
COPY --from=builder /opt/mega-xtea/cpp/extract_unmapped cpp/
COPY --from=builder /opt/mega-xtea/cpp/*.so cpp/
COPY mega-xtea.py .
COPY megaxtea/ megaxtea/
COPY scripts/ scripts/
COPY docs/ docs/

# Make entry point executable
RUN chmod +x mega-xtea.py && \
    ln -s /opt/mega-xtea/mega-xtea.py /usr/local/bin/mega-xtea

# Default data mount points
RUN mkdir -p /data/input /data/output /data/reference /data/kmer

# Verify installation
RUN python3 mega-xtea.py --version

ENTRYPOINT ["python3", "/opt/mega-xtea/mega-xtea.py"]
CMD ["--help"]

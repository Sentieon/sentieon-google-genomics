FROM google/cloud-sdk:333.0.0-slim as downloader

# Install samtools
RUN apt-get update && \
    apt-get install -y \
        lbzip2 \
        libncurses5-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        make && \
    curl -Lo samtools-1.9.tar.bz2 https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure && \
    make install

# Install samblaster
RUN curl -Lo samblaster-v.0.1.24.tar.gz https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.24/samblaster-v.0.1.24.tar.gz && \
  tar -xf samblaster-v.0.1.24.tar.gz && \
  cd samblaster-v.0.1.24 && \
  make && \
  cp samblaster /usr/local/bin/

FROM google/cloud-sdk:333.0.0-slim

LABEL container.base.image="google/cloud-sdk:333.0.0-slim"

COPY --from=downloader /usr/local/bin/samtools /usr/local/bin
COPY --from=downloader /usr/local/bin/samblaster /usr/local/bin

CMD ["/bin/bash"]

RUN apt-get update && apt-get install -y \
	curl \
        libncurses5-dev \
        bc \
        dnsutils \
        iputils-ping

# Install metadata script dependencies
RUN pip3 install requests urllib3

ADD gc_functions.sh gc_somatic.sh gc_germline.sh gc_ccdg_germline.sh gen_credentials.py /opt/sentieon/

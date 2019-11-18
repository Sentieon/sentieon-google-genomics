FROM python:2.7.16-slim-stretch as downloader

# Install gsutil
RUN apt-get update && \
    apt-get install -y \
        curl \
        python-pip \
        gcc \
        lsb-release && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    pip install crcmod

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

# Install metadata script dependencies
RUN pip install requests urllib3

FROM python:2.7.16-slim-stretch

LABEL container.base.image="python:2.7.16-slim-stretch"

COPY --from=downloader /usr/lib/google-cloud-sdk /usr/lib/google-cloud-sdk
COPY --from=downloader /usr/local/lib/python2.7/site-packages /usr/local/lib/python2.7/site-packages
COPY --from=downloader /usr/local/bin/samtools /usr/local/bin
COPY --from=downloader /usr/local/bin/samblaster /usr/local/bin

CMD ["/bin/bash"]

RUN apt-get update && apt-get install -y \
	curl \
        libncurses5-dev \
        bc \
        dnsutils \
        iputils-ping && \
    ln -s ../lib/google-cloud-sdk/bin/gsutil /usr/bin/gsutil

ADD gc_functions.sh gc_somatic.sh gc_germline.sh gc_ccdg_germline.sh gen_credentials.py /opt/sentieon/

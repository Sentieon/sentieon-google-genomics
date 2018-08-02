FROM sentieon/sentieon:201711.05 as builder
ARG version=201711.05

RUN apt-get update && \
    apt-get install -y \
        gcc \
        python-dev \
        python-setuptools \
        curl \
        lsb-release && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    pip install crcmod

RUN apt-get update && \
    apt-get install -y \
        lbzip2 \
        libncurses5-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        make && \
    curl -Lo samtools-1.7.tar.bz2 https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && \
    tar -xf samtools-1.7.tar.bz2 && \
    cd samtools-1.7 && \
    ./configure && \
    make install

FROM sentieon/sentieon:201711.05
ARG version=201711.05

COPY --from=builder /usr/lib/google-cloud-sdk /usr/lib/google-cloud-sdk
COPY --from=builder /usr/local/lib/python2.7/site-packages/crcmod/ /usr/local/lib/python2.7/site-packages/crcmod/
COPY --from=builder /usr/local/bin/samtools /usr/local/bin

RUN apt-get update && apt-get install -y \
        bc \
        dnsutils \
        iputils-ping && \
    pip install --upgrade requests && \
    ln -s ../lib/google-cloud-sdk/bin/gsutil /usr/bin/gsutil

ADD gc_functions.sh gc_somatic.sh gc_germline.sh gen_credentials.py poller.sh /opt/sentieon/


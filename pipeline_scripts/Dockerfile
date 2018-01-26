FROM sentieon/sentieon:201711.01 as gsutil-builder
ARG version=201711.01

RUN apt-get update && apt-get install -y \
    gcc \
    python-dev \
    python-setuptools \
    curl \
    lsb-release \
    && export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" \
    && echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - \
    && apt-get update && apt-get install -y google-cloud-sdk \
    && pip install crcmod

FROM sentieon/sentieon:201711.01
ARG version=201711.01

COPY --from=gsutil-builder /usr/lib/google-cloud-sdk /usr/lib/google-cloud-sdk
COPY --from=gsutil-builder /usr/local/lib/python2.7/site-packages/crcmod/ /usr/local/lib/python2.7/site-packages/crcmod/

RUN apt-get update && apt-get install -y \
    bc \
    && pip install --upgrade requests \
    && ln -s ../lib/google-cloud-sdk/bin/gsutil /usr/bin/gsutil

ADD gc_germline.sh gen_credentials.py /opt/sentieon/


FROM sentieon/sentieon:201711.01
ARG version=201711.01

RUN apt-get update && apt-get install -y \
    bc \
    lsb-release \
    curl \
    && export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" \
    && echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - \
    && apt-get update && apt-get install -y google-cloud-sdk \
    && pip install --upgrade requests

ADD gc_germline.sh gen_credentials.py /opt/sentieon/


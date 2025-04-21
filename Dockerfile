FROM amd64/python:3.9.9-slim-buster

USER root

# Install system dependencies
RUN apt-get update -y && \
	apt-get -y install build-essential wget unzip git

# Install PLINK v1.9 (equivalent to install-plink in Makefile)
ENV PLINK_DIR="/bin/plink"
RUN mkdir -p ${PLINK_DIR} && \
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip -O "${PLINK_DIR}/plink.zip" && \
    unzip "${PLINK_DIR}/plink.zip" -d "${PLINK_DIR}" && \
    rm "${PLINK_DIR}/plink.zip" && \
    chmod +x "${PLINK_DIR}/plink"
ENV PATH="${PLINK_DIR}:${PATH}"

# Install Python packages
RUN pip install snipar==0.0.20 dxpy
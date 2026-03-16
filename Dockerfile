#################################################################
# Dockerfile
#
# Version:          1.0
# Software:         OptiType
# Software Version: 1.3
# Description:      Accurate NGS-based 4-digit HLA typing
# Website:          https://github.com/FRED-2/OptiType/
# Tags:             Genomics
# Provides:         OptiType 1.3
# Base Image:       biodckr/biodocker
# Build Cmd:        docker build --rm -t fred2/opitype .
# Pull Cmd:         docker pull fred2/optitype
# Run Cmd:          docker run -v /path/to/file/dir:/data fred2/optitype
#################################################################

# Source Image
FROM ubuntu:16.04

################## BEGIN INSTALLATION ###########################
USER root

# install
RUN printf "deb http://dk.archive.ubuntu.com/ubuntu/ xenial main\ndeb http://dk.archive.ubuntu.com/ubuntu/ xenial universe\n" >> /etc/apt/sources.list
RUN apt-get update && apt-get install -y software-properties-common \
    gcc-4.9 \
    g++-4.9 \
    coinor-cbc \
    zlib1g-dev \
    libbz2-dev \
    libfreetype6-dev \
    libxft-dev \
    curl \
    git \
    cmake \
&& update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9 \
&& rm -rf /var/lib/apt/lists/* \
&& apt-get clean && apt-get purge

#HLA Typing
#OptiType dependecies
WORKDIR /
# COPY hdf5-1.8.21-Std-centos7-x86_64-shared_64.tar.gz .
RUN curl -O https://support.hdfgroup.org/ftp/HDF5/current18/bin/hdf5-1.8.21-Std-centos7-x86_64-shared_64.tar.gz && \
    tar -xvf hdf5-1.8.21-Std-centos7-x86_64-shared_64.tar.gz \
    && mv hdf5/bin/* /usr/local/bin/ \
    && mv hdf5/lib/* /usr/local/lib/ \
    && mv hdf5/include/* /usr/local/include/ \
    && mv hdf5/share/* /usr/local/share/ \
    && rm -rf hdf5/ \
    && rm -f hdf5-1.8.21-Std-centos7-x86_64-shared_64.tar.gz

ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH
ENV HDF5_DIR /usr/local/


RUN apt-get update && apt-get install -y python-setuptools python3-dev build-essential python3-pip && apt-get clean && apt-get purge

# pip install --upgrade pip &&
RUN pip3 install \
    numpy==1.10.0 \
    tables==3.2.2 \
    pysam==0.8.3 \
    future==0.15.2 \
    pyomo==4.2.10782 \
    pandas==0.23.0

RUN pip3 install matplotlib==1.4.3

# 虽然官方指出 Pandas 0.16.2,但安装报错 RuntimeError: Python version >= 3.5 required, 解决: 安装 pandas==0.23.0  ref:https://zhuanlan.zhihu.com/p/76850433
#     matplotlib==1.4.3 \
# 同理numpy==1.9.3 执行时报错 numpy==1.10.0

#installing optitype form git repository (version Dec 09 2015) and wirtig config.ini
RUN git clone https://github.com/FRED-2/OptiType.git \
    && sed -i -e '1i#!/usr/bin/env python' OptiType/OptiTypePipeline.py \
    && mv OptiType/ /usr/local/bin/ \
    && chmod 777 /usr/local/bin/OptiType/OptiTypePipeline.py \
    && echo "[mapping]\n\
razers3=/usr/local/bin/razers3 \n\
threads=1 \n\
\n\
[ilp]\n\
solver=cbc \n\
threads=1 \n\
\n\
[behavior]\n\
deletebam=true \n\
unpaired_weight=0 \n\
use_discordant=false\n" >> /usr/local/bin/OptiType/config.ini

#installing razers3
# git clone https://github.com/seqan/seqan.git seqan-src
WORKDIR /opt
RUN git clone https://github.com/seqan/seqan.git seqan-src
RUN cd seqan-src \
    && cmake -DCMAKE_BUILD_TYPE=Release \
    && make razers3 \
    && cp bin/razers3 /usr/local/bin/ \
    && cd .. \
    && rm -rf seqan-src

ENV PATH=/usr/local/bin/OptiType:$PATH


# Add user biodocker with password biodocker
RUN mkdir /data /config
RUN groupadd fuse && \
    useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo,fuse biodocker && \
    echo `echo "biodocker\nbiodocker\n" | passwd biodocker` && \
    chown biodocker:biodocker /data && \
    chown biodocker:biodocker /config

# Change user to back to biodocker
USER biodocker

# Change workdir to /data/
WORKDIR /data/

# Define default command
ENTRYPOINT ["python3", "OptiTypePipeline.py"]
CMD ["-h"]

##################### INSTALLATION END ##########################
# File Author / Maintainer
MAINTAINER Benjamin Schubert <schubert@informatik.uni-tuebingen.de>
# modify by biolxy <biolxy@aliyun.com>
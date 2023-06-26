FROM ubuntu:18.04
LABEL Chosen Obih
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y g++ \
		make \
		git \
		zlib1g-dev \
		openssl \
		lbzip2 \
		bzip2 \
		python \
		python2.7-dev \
		perl \
		wget \
		curl \
		python-matplotlib \
		python-numpy \
        python-pandas

ENV BINPATH /usr/bin
WORKDIR /evolinc_docker

# Cufflinks
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -

# Transdecoder
RUN wget -O- https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz | tar xzvf -

# Diamond Blast
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.10/diamond-linux64.tar.gz
RUN tar xzf diamond-linux64.tar.gz

# Samtools
RUN wget --no-check-certificate http://sourceforge.net/projects/samtools/files/samtools/1.0/samtools-bcftools-htslib-1.0_x64-linux.tar.bz2/download
RUN tar xvf download

# Bedtools
RUN wget https://github.com/arq5x/bedtools2/archive/v2.25.0.tar.gz
RUN tar xvf v2.25.0.tar.gz
RUN cd bedtools2-2.25.0 && make
RUN cd ..

# Bedops tool
RUN wget -O- https://github.com/bedops/bedops/releases/download/v2.4.16/bedops_linux_x86_64-v2.4.16.tar.bz2 | tar jxvf -

# cpan
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm URI/Escape.pm

# R libraries
RUN apt-get -y install ca-certificates software-properties-common gnupg2 gnupg1 gnupg
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/"
RUN apt-get install -y r-base
RUN Rscript -e 'install.packages("splitstackshape", dependencies = TRUE, repos="http://cran.rstudio.com/");'
RUN Rscript -e 'install.packages("dplyr", dependencies = TRUE, repos="http://cran.rstudio.com/");'
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")};'
RUN Rscript -e 'BiocManager::install(c("Biostrings"));'
RUN Rscript -e 'install.packages("getopt", dependencies = TRUE, repos="http://cran.rstudio.com/");'

# Uniprot database
ADD https://github.com/iPlantCollaborativeOpenSource/docker-builds/releases/download/evolinc-I/uniprot_sprot.dmnd.gz /evolinc_docker/
RUN gzip -d /evolinc_docker/uniprot_sprot.dmnd.gz

# rFAM database
ADD https://de.cyverse.org/dl/d/12EF1A2F-B9FC-456D-8CD9-9F87197CACF2/rFAM_sequences.fasta /evolinc_docker/

# Biopython
RUN curl "https://bootstrap.pypa.io/pip/2.7/get-pip.py" -o "get-pip.py"
RUN curl "https://files.pythonhosted.org/packages/27/79/8a850fe3496446ff0d584327ae44e7500daf6764ca1a382d2d02789accf7/pip-20.3.4-py2.py3-none-any.whl" -o "pip-20.3.4-py2.py3-none-any.whl"
RUN curl "https://files.pythonhosted.org/packages/e1/b7/182161210a13158cd3ccc41ee19aadef54496b74f2817cc147006ec932b4/setuptools-44.1.1-py2.py3-none-any.whl" -o "setuptools-44.1.1-py2.py3-none-any.whl"
RUN curl "https://files.pythonhosted.org/packages/27/d6/003e593296a85fd6ed616ed962795b2f87709c3eee2bca4f6d0fe55c6d00/wheel-0.37.1-py2.py3-none-any.whl" -o "wheel-0.37.1-py2.py3-none-any.whl"
RUN python get-pip.py "pip-20.3.4-py2.py3-none-any.whl" "setuptools-44.1.1-py2.py3-none-any.whl" "wheel-0.37.1-py2.py3-none-any.whl"
RUN curl "https://files.pythonhosted.org/packages/ff/f4/0ce39bebcbb0ff619426f2bbe86e60bc549ace318c5a9113ae480ab2adc7/biopython-1.76.tar.gz" -o "biopython-1.76.tar.gz"
RUN tar xvf biopython-1.76.tar.gz
RUN cd biopython-1.76 && python setup.py build && python setup.py install
RUN cd ..

# CPC2
ADD CPC2-beta /evolinc_docker/CPC2-beta
WORKDIR /evolinc_docker/CPC2-beta/libs/libsvm/
RUN tar xvf libsvm-3.22.tar.gz
WORKDIR libsvm-3.22
RUN make clean && make
WORKDIR /

#LAST
RUN wget https://gitlab.com/mcfrith/last/-/archive/main/last-main.zip
RUN apt-get install -y unzip
RUN unzip last-main.zip
WORKDIR last-main
RUN make
RUN make install prefix=~

# Evolinc wrapper scripts
ADD *.sh *.py *.R /evolinc_docker/
RUN chmod +x /evolinc_docker/evolinc-part-I.sh && cp /evolinc_docker/evolinc-part-I.sh $BINPATH
WORKDIR /

# Setting paths to all the softwares
ENV PATH /evolinc_docker/cufflinks-2.2.1.Linux_x86_64/:$PATH
ENV PATH /evolinc_docker/TransDecoder-2.0.1/:$PATH
ENV PATH /evolinc_docker/ncbi-blast-2.4.0+/bin/:$PATH
ENV PATH /evolinc_docker/bedtools2-2.25.0/bin/:$PATH
ENV PATH /evolinc_docker/samtools-bcftools-htslib-1.0_x64-linux/bin/:$PATH
ENV PATH /evolinc_docker/bin/:$PATH
ENV PATH /evolinc_docker/CPC2-beta/bin/:$PATH
ENV PATH /last-main/src/:$PATH
ENV PATH /evolinc_docker/:$PATH

# Entrypoint
ENTRYPOINT ["evolinc-part-I.sh"]
CMD ["-h"]
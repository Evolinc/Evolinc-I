FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty

RUN apt-get update && apt-get install -y g++ \
		make \
		git \
		zlib1g-dev \
		python \
		wget \
		curl \
		python-matplotlib \
		python-numpy \
                python-pandas

ENV BINPATH /usr/bin
ENV EVGIT https://upendra_35@bitbucket.org/upendra_35/evolinc_docker.git

RUN git clone $EVGIT
WORKDIR /evolinc_docker

RUN chmod +x evolinc-part-I.sh && cp evolinc-part-I.sh $BINPATH
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -
RUN wget -O- https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz | tar xzvf -
RUN wget -O- http://seq.cs.iastate.edu/CAP3/cap3.linux.x86_64.tar | tar vfx -
RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.3.0+-x64-linux.tar.gz > ncbi-blast-2.3.0+-x64-linux.tar.gz
RUN tar xvf ncbi-blast-2.3.0+-x64-linux.tar.gz
RUN wget -O- http://ftp.mirrorservice.org/sites/download.sourceforge.net/pub/sourceforge/q/qu/quast/quast-3.0.tar.gz | tar zxvf -
RUN wget http://sourceforge.net/projects/samtools/files/samtools/1.0/samtools-bcftools-htslib-1.0_x64-linux.tar.bz2/download
RUN tar xvf download
RUN wget https://github.com/arq5x/bedtools2/archive/v2.25.0.tar.gz
RUN tar xvf v2.25.0.tar.gz
RUN cd bedtools2-2.25.0 && make
RUN cd ..
RUN wget http://sourceforge.net/projects/bio-bwa/files/latest/download?source=files
RUN tar xvf download?source=files
RUN cd bwa-0.7.12 && make
RUN cd ..
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm URI/Escape.pm

ENV PATH /evolinc_docker/cufflinks-2.2.1.Linux_x86_64/:$PATH
ENV PATH /evolinc_docker/TransDecoder-2.0.1/:$PATH
ENV PATH /evolinc_docker/CAP3/:$PATH
ENV PATH /evolinc_docker/ncbi-blast-2.3.0+/bin/:$PATH
ENV PATH /evolinc_docker/bedtools2-2.25.0/bin/:$PATH
ENV PATH /evolinc_docker/samtools-bcftools-htslib-1.0_x64-linux/bin/:$PATH
ENV PATH /evolinc_docker/bwa-0.7.12/:$PATH

ENTRYPOINT ["/usr/bin/evolinc-part-I.sh"]
CMD ["-h"]

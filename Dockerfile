FROM ubuntu:14.04.3
MAINTAINER Upendra Devisetty

RUN apt-get update && apt-get install -y g++ \
		make \
		git \
		zlib1g-dev \
		python \
		wget \
		curl \
		python-matplotlib

ENV BINPATH /usr/bin
ENV EVO2GIT https://upendra_35@bitbucket.org/upendra_35/evolinc_docker.git

RUN git clone $EVO2GIT
RUN chmod +x /evolinc_docker/evolinc-part-I.sh && cp /evolinc_docker/evolinc-part-I.sh $BINPATH
RUN wget -O- http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz | tar xzvf -
RUN wget -O- https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz | tar xzvf -
RUN wget -O- http://seq.cs.iastate.edu/CAP3/cap3.linux.x86_64.tar | tar vfx -
RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.31+-x64-linux.tar.gz > ncbi-blast-2.2.31+-x64-linux.tar.gz
RUN tar xvf ncbi-blast-2.2.31+-x64-linux.tar.gz
RUN wget -O- http://ftp.mirrorservice.org/sites/download.sourceforge.net/pub/sourceforge/q/qu/quast/quast-3.0.tar.gz | tar zxvf -
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm URI/Escape.pm

ENV PATH /CAP3/:$PATH
ENV PATH /ncbi-blast-2.2.31+/bin/:$PATH
ENV PATH /quast-3.0/:$PATH
ENV PATH /cufflinks-2.2.1.Linux_x86_64/:$PATH
ENV PATH /TransDecoder-2.0.1/:$PATH

ENTRYPOINT ["/usr/bin/evolinc-part-I.sh"]
CMD ["-h"]

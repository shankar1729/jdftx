FROM ubuntu:16.04
RUN apt-get -y update && apt-get -y --no-install-recommends install g++ cmake libgsl0-dev libopenmpi-dev openmpi-bin libfftw3-dev libatlas-base-dev liblapack-dev wget unzip ca-certificates make && \
	cd /root && \
	wget https://github.com/shankar1729/jdftx/archive/refs/heads/master.zip && unzip master.zip && rm master.zip && \
	cd jdftx-master && mkdir build && cd build && \
	cmake ../jdftx && make all && make install && \
        make test && \
        apt-get -y purge gcc g++ cmake cmake-data g++-5 gfortran-5 make libarchive13 libcurl3 libgfortran-5-dev libjsoncpp1 liblzo2-2 libstdc++-5-dev wget unzip ca-certificates && rm -rf /var/lib/apt/lists/* && \
        cd /root && rm -rf /root/jdftx-master && \
        echo 'export PATH="$PATH:/usr/local/share/jdftx/scripts"' >> /root/.bashrc && mkdir /root/research
WORKDIR /root/research

#Use it like this:
#docker run -it --rm -v $PWD:/root/research jdftx

#Or even better, put the following line at the end of your .bashrc and/or .zshrc so that you can just run 'jdftx' :
#function jdftx () { docker run -it --rm -v $PWD:/root/research jdftx ; }

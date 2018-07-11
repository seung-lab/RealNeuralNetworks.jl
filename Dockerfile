FROM julia:0.6
LABEL   maintainer="Jingpeng Wu <jingpeng@princeton.edu>" \
        function="skeletonization"

RUN apt-get update
RUN apt-get install -qq --no-install-recommends build-essential unzip wget \
                libjemalloc-dev libhdf5-dev libxml2 libmagickwand-dev
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so 

RUN julia -e 'Pkg.update()'
RUN julia -e 'Pkg.clone("https://github.com/jingpengw/HTTP.jl.git")'
RUN julia -e 'Pkg.add("LightGraphs")'
RUN julia -e 'Pkg.add("MetaGraphs")'
RUN julia -e 'Pkg.add("ArgParse")'
RUN julia -e 'Pkg.add("LsqFit")'
RUN julia -e 'Pkg.add("DataFrames")'
RUN julia -e 'Pkg.add("ImageFiltering")'
RUN julia -e 'Pkg.add("OffsetArrays")'
RUN julia -e 'Pkg.add("Query")'
# FFTW was needed to fix a julia version problem
# https://github.com/JuliaMath/FFTW.jl/issues/45
RUN julia -e 'Pkg.add("FFTW")'
RUN julia -e 'Pkg.clone("https://github.com/jingpengw/AWSCore.jl.git")'
RUN julia -e 'Pkg.clone("https://github.com/jingpengw/AWSSDK.jl.git")'
RUN julia -e 'Pkg.clone("https://github.com/jingpengw/AWSS3.jl.git")'
RUN julia -e 'Pkg.clone("https://github.com/seung-lab/EMIRT.jl.git")'
RUN julia -e 'Pkg.clone("https://github.com/seung-lab/BigArrays.jl.git")'

WORKDIR /root/.julia/v0.6/
RUN mkdir RealNeuralNetworks 
ADD . RealNeuralNetworks/ 

RUN julia -e 'using RealNeuralNetworks'
WORKDIR /root/.julia/v0.6/RealNeuralNetworks/scripts
#CMD ["julia", "skeletonize.jl", "-h"]

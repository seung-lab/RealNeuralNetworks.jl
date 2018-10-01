FROM julia:0.7
LABEL   maintainer="Jingpeng Wu <jingpeng@princeton.edu>" \
        function="skeletonization"

RUN apt-get update
RUN apt-get install -qq --no-install-recommends build-essential unzip wget \
                libjemalloc-dev libhdf5-dev libxml2 libmagickwand-dev
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so 

RUN julia -e 'import Pkg; Pkg.update()'

RUN julia -e 'import Pkg; \ 
                Pkg.add("LightGraphs");  \
                Pkg.add("MetaGraphs"); \
                Pkg.add("ArgParse");\
                Pkg.add("LsqFit");\
                Pkg.add("DataFrames");\
                Pkg.add("ImageFiltering");\
                Pkg.add("OffsetArrays"); \
                Pkg.add("Query")'

#RUN julia -e 'import Pkg; Pkg.clone("https://github.com/seung-lab/EMIRT.jl.git")'
RUN julia -e 'import Pkg; Pkg.clone("https://github.com/seung-lab/BigArrays.jl.git")'
RUN julia -e 'import Pkg; Pkg.clone("https://github.com/seung-lab/RealNeuralNetworks.jl.git")'


#WORKDIR /root/.julia/dev/
#RUN mkdir RealNeuralNetworks 
#ADD . RealNeuralNetworks/ 
#RUN julia -e 'import Pkg; Pkg.clone(pwd()); Pkg.test("RealNeuralNetworks")'

RUN julia -e 'using RealNeuralNetworks'
WORKDIR /root/.julia/packages/RealNeuralNetworks/scripts
#CMD ["julia", "skeletonize.jl", "-h"]

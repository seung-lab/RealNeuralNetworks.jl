FROM julia:1.3
LABEL   maintainer="Jingpeng Wu <jingpeng@princeton.edu>" \
        function="skeletonization"

RUN apt-get update
RUN apt-get install -qq --no-install-recommends build-essential unzip wget \
                libjemalloc-dev libhdf5-dev libxml2 libmagickwand-dev
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so 

WORKDIR /root/.julia/dev/
RUN mkdir RealNeuralNetworks 
ADD . RealNeuralNetworks/ 


# install latest packages first to make sure that we are using the latest 
# RUN julia -e 'import Pkg; Pkg.clone("https://github.com/seung-lab/EMIRT.jl.git")'
#RUN julia -e 'import Pkg; Pkg.clone("https://github.com/JuliaCloud/GoogleCloud.jl.git")'
RUN julia -e 'using Pkg; Pkg.update(); Pkg.instantiate();' 
    #Pkg.add("ImageMagick"); using ImageMagick; \
RUN julia -e 'using Pkg; Pkg.develop(PackageSpec(name="BigArrays", url="https://github.com/seung-lab/BigArrays.jl.git"))'
    #Pkg.resolve(); \
    #Pkg.develop(PackageSpec(url="https://github.com/seung-lab/RealNeuralNetworks.jl.git")); \
RUN julia -e 'using Pkg; Pkg.develop(PackageSpec(name="RealNeuralNetworks", path=pwd())); \
        Pkg.develop("RealNeuralNetworks"); Pkg.instantiate();' 
    # install registered packages later
RUN julia -e 'using Pkg; Pkg.add("LightGraphs");  \
    Pkg.add("MetaGraphs"); \
    Pkg.add("ArgParse");\
    Pkg.add("LsqFit");\
    Pkg.add("DataFrames");\
    Pkg.add("ImageFiltering");\
    Pkg.add("OffsetArrays"); \
    Pkg.add("Query"); \
    Pkg.add("JSON");\ 
    Pkg.add("AWSSDK"); '

# https://discourse.julialang.org/t/pkg-add-ijulia-can-not-work/13341/2
#RUN julia -e 'rm(joinpath(homedir(), ".julia", "registries"); recursive=true)'

# precompile the package
RUN julia -e 'using RealNeuralNetworks'
WORKDIR /root/.julia/dev/RealNeuralNetworks/scripts
#CMD ["julia", "skeletonize.jl", "-h"]

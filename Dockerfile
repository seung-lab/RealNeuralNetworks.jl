FROM julia:1.0
LABEL   maintainer="Jingpeng Wu <jingpeng@princeton.edu>" \
        function="skeletonization"

RUN apt-get update
RUN apt-get install -qq --no-install-recommends build-essential unzip wget \
                libjemalloc-dev libhdf5-dev libxml2 libmagickwand-dev
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so 

RUN julia -e 'import Pkg; Pkg.update()'

# install latest packages first to make sure that we are using the latest 
# RUN julia -e 'import Pkg; Pkg.clone("https://github.com/seung-lab/EMIRT.jl.git")'
#RUN julia -e 'import Pkg; Pkg.clone("https://github.com/JuliaCloud/GoogleCloud.jl.git")'
RUN julia -e 'using Pkg; Pkg.instantiate();\
    #Pkg.add("ImageMagick"); using ImageMagick; \
    Pkg.develop(PackageSpec(name="BigArrays", url="https://github.com/seung-lab/BigArrays.jl.git")); \
    Pkg.resolve(); \
    #Pkg.develop(PackageSpec(url="https://github.com/seung-lab/RealNeuralNetworks.jl.git")); \
    Pkg.develop(PackageSpec(name="RealNeuralNetworks", path=pwd())); \
    # install registered packages later
    Pkg.add("LightGraphs");  \
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

#WORKDIR /root/.julia/dev/
#RUN mkdir RealNeuralNetworks 
#ADD . RealNeuralNetworks/ 
#RUN julia -e 'import Pkg; Pkg.clone(pwd()); Pkg.test("RealNeuralNetworks")'

# precompile the package
RUN julia -e 'using RealNeuralNetworks'
WORKDIR /root/.julia/dev/RealNeuralNetworks/scripts
#CMD ["julia", "skeletonize.jl", "-h"]

RealNeuralNetworks.jl
========================
3D neuron models extracted from EM image segmentation

## Installation
run this inside julia REPL:

    Pkg.clone("https://github.com/seung-lab/RealNeuralNetworks.jl.git")

## Usage
you must have google secret json file in `/secrets/google-secret.json`, same with configuration of [GSDicts.jl](https://github.com/seung-lab/GSDicts.jl) or [cloudvolume](https://github.com/seung-lab/cloud-volume)

### commandline
`julia ~/julia/v0.5/RealNeuralNetworks/scripts/skeletonize.jl -h`

### Docker (recommanded)
#### build docker image
    cd ~/julia/v0.5/RealNeuralNetworks
    sudo docker build . -t realneuralnetworks

```
docker run -v /tmp:/tmp -v /secrets:/secrets --net=host realneuralnetworks julia skeletonize.jl -h
```

### REPL in Julia

```Julia
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.BranchNets
using RealNeuralNetworks.SWCs

nodeNet = NodeNet(seg::Array{UInt32,3}; obj_id = convert(UInt32,77605))
branchNet = BranchNet( nodeNet )
swc = SWC(branchNet)
SWCs.save(swc, tempname()*".swc")
```

## Credit 
The skeletonization was originally implemented in Matlab by Alexander Bae using TEASAR algorithm, which was translated to Julia by Nicholas Turner.

## Skeletonization Algorithm 
Sato, M., I. Bitter, M. A. Bender, A. E. Kaufman, and M. Nakajima. “TEASAR: Tree-Structure Extraction Algorithm for Accurate and Robust Skeletons.” In Proceedings the Eighth Pacific Conference on Computer Graphics and Applications, 281–449, 2000. doi:10.1109/PCCGA.2000.883951.



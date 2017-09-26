TEASAR.jl
========================
3D chunk-wise skeletonization using TEASAR algorithm

## Installation
run this inside julia REPL:

    Pkg.clone("https://github.com/seung-lab/TEASAR.jl.git")

## Usage
you must have google secret json file in `/secrets/google-secret.json`, same with configuration of [GSDicts.jl](https://github.com/seung-lab/GSDicts.jl) or [cloudvolume](https://github.com/seung-lab/cloud-volume)

### commandline
`julia ~/julia/v0.5/TEASAR/scripts/skeletonize.jl -h`

### Docker (recommanded)
#### build docker image
    cd ~/julia/v0.5/TEASAR
    sudo docker build . -t teasar

```
docker run -v /tmp:/tmp -v /secrets:/secrets --net=host teasar julia skeletonize.jl -h
```

### REPL in Julia

```Julia
using TEASAR.Skeletons
using TEASAR.SWCs

skeleton = Skeleton(seg::Array{UInt32,3}; obj_id = convert(UInt32,77605))
swc = SWC(skeleton)
TEASAR.SWCs.save(swc, tempname()*".swc")
```
or
```
skeleton = Skeleton(points::Array{Int,2})
```

## Algorithm 
Sato, M., I. Bitter, M. A. Bender, A. E. Kaufman, and M. Nakajima. “TEASAR: Tree-Structure Extraction Algorithm for Accurate and Robust Skeletons.” In Proceedings the Eighth Pacific Conference on Computer Graphics and Applications, 281–449, 2000. doi:10.1109/PCCGA.2000.883951.

## Credit 
originally implemented by Alexander Bae using matlab, then translated to Julia by Nicholas Turner.

# TEASAR.jl
3D chunk-wise skeletonization using TEASAR algorithm

# Installation
run this inside julia REPL:

    Pkg.clone("https://github.com/seung-lab/TEASAR.jl.git")

# Usage
```Julia
using TEASAR
Teasar(points)
```

# Algorithm 
Sato, M., I. Bitter, M. A. Bender, A. E. Kaufman, and M. Nakajima. “TEASAR: Tree-Structure Extraction Algorithm for Accurate and Robust Skeletons.” In Proceedings the Eighth Pacific Conference on Computer Graphics and Applications, 281–449, 2000. doi:10.1109/PCCGA.2000.883951.

# Credit 
originally implemented by Alexander Bae using matlab, then translated to Julia by Nicholas Turner.

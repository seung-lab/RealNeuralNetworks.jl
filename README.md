RealNeuralNetworks.jl
========================

[![Build Status](https://travis-ci.org/seung-lab/RealNeuralNetworks.jl.svg?branch=master)](https://travis-ci.org/seung-lab/RealNeuralNetworks.jl)

3D neuron models extracted from EM image segmentation

# Installation
run this inside julia package REPL (hit `]` to enter package mode):

    dev https://github.com/seung-lab/RealNeuralNetworks.jl.git 

# Usage

## commandline
`julia skeletonize.jl -h`

## Docker (recommanded)
### build docker image
    sudo docker build . -t realneuralnetworks

```
docker run -v /tmp:/tmp -v /secrets:/secrets --net=host realneuralnetworks julia skeletonize.jl -h
```

## REPL in Julia

```Julia
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.Neurons
using RealNeuralNetworks.SWCs

nodeNet = NodeNet(seg::Array{UInt32,3}; obj_id = convert(UInt32,77605))
neuron = Neuron( nodeNet )
swc = SWC(neuron)
SWCs.save(swc, tempname()*".swc")
```

# Morphological Features

## features represent a whole neuron

- [x] arbor density
- [x] total path length 
- [x] number of segment points 
- [x] Median segment length is the median dendritic segment length of all the segments starting and ending at irreducible nodes (in μm). Irreducible nodes are the points of the dendritic arbor corresponding to soma, branching points or terminal points.
- [x] 3D sholl analysis. 
- [ ] Hull area is the area of the tightest convex hull containing the z-projection of the dendritic arbor (in μm2). 
- [ ] volume of the convex hull around all neurites
- [x] Average angle is the mean of the positive angle between (parent node, node) and (node, child node) vectors, where node, parent node and child node are irreducible nodes (in radians).  
- [x] Average tortuosity is the average value of the ratio of the actual dendritic length to the Euclidean distance between irreducible nodes. 
- [x] Asymmetry is the distance of the soma node to the dendritic arbor (skeleton) centre of mass (in nm). 
- [x] Typical radius (λ) is the root-mean-square distance of dendritic arbor points to the centre of mass (in nm). 
- [x] fractal dimension.
- [x] longest neurite length.
- [ ] [Strahler number](https://en.wikipedia.org/wiki/Strahler_number)
- [ ] distribution of Euclidian distance of segmentes from soma (third principal component)
- [ ] distribution of Euclidian distance of segmentes from soma as a function of segment order (third principal component)
- [ ] number of segmentes per segment order (second principal component)
- [ ] distribution of morphological distance of segmentes from soma along the skeleton as a function of segment order (first principal component)

## features represent segmentes in a neuron
- [x] ratio of tail diameter to head. could be useful to identify spines. 
- [x] segment order
- [x] segment length
- [x] branching angle. [computation using dot product](https://stackoverflow.com/questions/19729831/angle-between-3-points-in-3d-space)
- [x] tortuosity / curvature. caa4486b501f936743f781782e6561833da7e413
- [x] distance to root path length
- [ ] [segment asymmetry](http://www.treestoolbox.org/manual/asym_tree.html)
- [x] average radius. easy to compute with radius list.

# Credit 
The skeletonization was originally implemented in Matlab by Alexander Bae using TEASAR algorithm, which was translated to Julia by Nicholas Turner.

# References:
1. Sümbül U, Song S, McCulloch K, Becker M, Lin B, Sanes JR, Masland RH, Seung HS. A genetic and computational approach to structurally classify neuronal types. Nature communications. 2014 Mar 24;5:3512. [link](https://www.nature.com/articles/ncomms4512#methods)
1. Schierwagen A, Villmann T, Alpár A, Gärtner U. Cluster analysis of cortical pyramidal neurons using som. InIAPR Workshop on Artificial Neural Networks in Pattern Recognition 2010 Apr 11 (pp. 120-130). Springer, Berlin, Heidelberg.
1. Cuntz H, Forstner F, Borst A, H\äusser M. The TREES Toolbox—Probing the Basis of Axonal and Dendritic Segmenting. Neuroinformatics. 2011;1–6. 
1. Schierwagen A. Neuronal morphology: Shape characteristics and models. Neurophysiology. 2008;40(4):310–315. 
1. Uylings HB., van Pelt J. Measures for quantifying dendritic arborizations. Network: Computation in Neural Systems. 2002;13(3):397–414. [a good review paper](http://www.tandfonline.com/doi/abs/10.1088/0954-898X_13_3_309)
1. Wanner AA, Genoud C, Masudi T, Siksou L, Friedrich RW. Dense EM-based reconstruction of the interglomerular projectome in the zebrafish olfactory bulb. Nature neuroscience. 2016 Jun 1;19(6):816-25. [also have clustering methods](https://www.nature.com/neuro/journal/v19/n6/abs/nn.4290.html)
1. Sato, M., I. Bitter, M. A. Bender, A. E. Kaufman, and M. Nakajima. “TEASAR: Tree-Structure Extraction Algorithm for Accurate and Robust Skeletons.” In Proceedings the Eighth Pacific Conference on Computer Graphics and Applications, 281–449, 2000. doi:10.1109/PCCGA.2000.883951.

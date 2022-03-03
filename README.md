![](https://github.com/seung-lab/RealNeuralNetworks.jl/blob/master/asset/logo/RGB_for_web/RNN_Logo_RGB.jpg)
========================
[![Generic badge](https://img.shields.io/badge/Documentation-passing.svg)](https://seung-lab.github.io/RealNeuralNetworks.jl/latest/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://seung-lab.github.io/RealNeuralNetworks.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://seung-lab.github.io/RealNeuralNetworks.jl/dev)
[![CI](https://github.com/seung-lab/RealNeuralNetworks.jl/workflows/CI/badge.svg)](https://github.com/seung-lab/RealNeuralNetworks.jl/actions?query=workflow%3ACI)
[![Coverage Status](https://coveralls.io/repos/github/seung-lab/RealNeuralNetworks.jl/badge.svg?branch=master)](https://coveralls.io/github/seung-lab/RealNeuralNetworks.jl?branch=master)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
<!-- [![Build Status](https://travis-ci.org/seung-lab/RealNeuralNetworks.jl.svg?branch=master)](https://travis-ci.org/seung-lab/RealNeuralNetworks.jl) -->


3D neuron models extracted from EM image segmentation 

# Features 
- skeletonization. extract the neuron skeletons based on image segmentation (colored labels).
- morphological analysis with a lot of features. 
- neural networks including synapses. Most of current morphological analysis tools do not have synapses. 
- cell type classification based on NBLAST. use NBLAST to compute similarity scores between neurons. 

# Installation
This package is registered in Julia package management system. Simply following these steps to install it.
- type `julia` to enter command line interface.
- inside julia package REPL, hit `]` to enter package mode, then type `add RealNeuralNetworks`

# Usage

## Analysis
All the analysis functions are demonstrated in [Jupyter Notebooks](https://github.com/jingpengw/realneuralnetworks-notebook).

## skeletonization of proofread neurons
The skeletonization could be done using a script. The following command will show the parameters to use it.
`julia skeletonize.jl -h`

## Docker
### build docker image
    sudo docker build . -t realneuralnetworks

The secret files containing cloud authentication keys follows the [cloudvolume secret format](https://github.com/seung-lab/cloud-volume#credentials).

```
docker run -v /tmp:/tmp -v $HOME/.cloudvolume/secrets:/root/.cloudvolume/secrets --net=host realneuralnetworks julia skeletonize.jl -h
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
- [x] Hull area is the area of the tightest convex hull containing the z-projection of the dendritic arbor (in μm2). 
- [x] volume of the convex hull around all neurites
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

# Reference
We have a [publication](https://www.frontiersin.org/articles/10.3389/fninf.2022.828169/full) for this repo:
```bibtex

@article{wuRealNeuralNetworksJlIntegrated2022,
	title = {{RealNeuralNetworks}.jl: {An} {Integrated} {Julia} {Package} for {Skeletonization}, {Morphological} {Analysis}, and {Synaptic} {Connectivity} {Analysis} of {Terabyte}-{Scale} {3D} {Neural} {Segmentations}},
	volume = {16},
	issn = {1662-5196},
	shorttitle = {{RealNeuralNetworks}.jl},
	url = {https://www.frontiersin.org/article/10.3389/fninf.2022.828169},
	urldate = {2022-03-02},
	journal = {Frontiers in Neuroinformatics},
	author = {Wu, Jingpeng and Turner, Nicholas and Bae, J. Alexander and Vishwanathan, Ashwin and Seung, H. Sebastian},
	year = {2022},
}

```

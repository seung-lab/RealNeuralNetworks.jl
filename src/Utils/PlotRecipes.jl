module PlotRecipes
using RealNeuralNetworks.Neurons
using Colors, ColorSchemes, Clustering
using Plots
using PyPlot

@everywhere const VOXEL_SIZE = (400,400,400)
@everywhere const GAUSSIAN_FILTER_STD = 8.0

function reorder_distance_matrix(distanceMatrix::Array, clust::Hclust)
	# rearrange distance matrix according to order
	reorderedDistanceMatrix = zeros(distanceMatrix)
	for (i,o1) in enumerate(clust.order)
		for (j, o2) in enumerate(clust.order)
			reorderedDistanceMatrix[i,j] = distanceMatrix[o1,o2]
		end 
	end 
	reorderedDistanceMatrix
end 

function coloring(dm)
    img = Array{RGB{Float64}}(size(dm))
    for i in eachindex(dm)
        img[i] = get(ColorSchemes.jet, dm[i])
    end 
    img
end 

function plot(neuron::Neuron; nodeStep = 10)
    segmentList = neuron.segmentList
    plotly()
    for branch in segmentList
        nodeList = Neurons.Segments.get_node_list(branch)
        #@show length(nodeList)
        #@show nodeList
        x = map(n->n[1], nodeList[1:nodeStep:end])
        y = map(n->n[2], nodeList[1:nodeStep:end])
        z = map(n->n[3], nodeList[1:nodeStep:end])
        plot!(x,y,z) #color=rand(Colors.RGB))
    end 

    root = Neurons.get_root_node(neuron)
    plot!([root[1]], [root[2]], [root[3]], m=(2, :circle), leg=false)
end 

function plot_arbor_density_map(densityMap::Array)
    fig = PyPlot.figure()
    fig[:add_subplot](2,2,1)
    xy = maximum(densityMap, 3)[:,:,1]
    PyPlot.imshow(xy, "jet")
    # colorbar()
    fig[:add_subplot](2,2,2)
    xz = maximum(densityMap, 2)[:,1,:]
    PyPlot.imshow(xz, "jet")
    fig[:add_subplot](2,2,3)
    yz = maximum(densityMap, 1)[1,:,:] |> rotl90
    PyPlot.imshow(yz, "jet")
    colorbar()
end 

function plot_mask(mask::Array{T,3}) where T
    fig = PyPlot.figure()
    fig[:add_subplot](2,2,1)
    xy = maximum(mask, 3)[:,:,1]
    PyPlot.imshow(xy, "gray")
    # colorbar()
    fig[:add_subplot](2,2,2)
    xz = maximum(mask, 2)[:,1,:]
    PyPlot.imshow(xz, "gray")
    fig[:add_subplot](2,2,3)
    yz = maximum(mask, 1)[1,:,:] |> rotl90
    PyPlot.imshow(yz, "gray")
end 
 
function indexmap(x::Vector)
    ret = Dict()
    for (i,v) in enumerate(x)
        ret[v] = i
    end 
    ret
end
function hclustplot(hc::Hclust, useheight::Bool)
    o = indexmap(hc.order)
    n = [x for x in 1:length(o)]

    pos = treepositions(hc, useheight)


    xs = []
    ys = []
    for i in 1: size(hc.merge, 1)
        x1 = pos[hc.merge[i,1]][1]
        x2 = pos[hc.merge[i,2]][1]
        append!(xs, [x1,x1,x2,x2])

        y1 = pos[hc.merge[i,1]][2]
        y2 = pos[hc.merge[i,2]][2]
        useheight ? h = hc.height[i] : h = 1
        newy = maximum([y1,y2]) + h
        append!(ys, [y1,newy,newy,y2])
    end
    return (reshape(xs, 4, size(hc.merge, 1)), reshape(ys, 4, size(hc.merge, 1)))
end

function treepositions(hc::Hclust, useheight::Bool)
    order = indexmap(hc.order)
    positions = Dict{}()
    for (k,v) in order
        positions[-k] = (v, 0)
    end
    for i in 1:size(hc.merge,1)
        xpos = mean([positions[hc.merge[i,1]][1], positions[hc.merge[i,2]][1]])
        if hc.merge[i,1] < 0 && hc.merge[i,2] < 0
            useheight ? ypos = hc.height[i] : ypos = 1
        else
            useheight ? h = hc.height[i] : h = 1
            ypos = maximum([positions[hc.merge[i,1]][2], positions[hc.merge[i,2]][2]]) + h
        end

        positions[i] = (xpos, ypos)
    end
    return positions
end 

function plot(clust::Hclust)
	Plots.plotly()
	Plots.plot(hclustplot(clust, true), seriestype=:path, color=:black,
    yaxis=nothing,  grid=false, legend=false) #,  xticks=classificationIdList[clust.order])
end 

end # module

module PlotRecipes

import PyPlot
PyPlot.svg(true)
using Plots
using Statistics 
#using Plotly
using PlotlyJS
using Clustering
using SparseArrays

using RealNeuralNetworks.Neurons
using RealNeuralNetworks.Neurons.Segments
using RealNeuralNetworks.Neurons.Segments.Synapses

function plot_synapse_distributions( cellList;  
        synapseDistribution = map(get_synapse_to_soma_path_length_lists, cellList),
        saveName = joinpath(homedir(), "synapse_distribution.svg"),
        prerange=(0,50), postrange=(0,100), bins=30)
    fig = PyPlot.figure()
    for (i,cellId) in enumerate(cellList)
        @show i
        fig[:add_subplot](4,2,(i-1)*2+1)
        PyPlot.plt[:hist](synapseDistribution[i][1]./1000; range=prerange, bins=bins)
        PyPlot.plt[:xlabel]("distrance from presynapse to soma ($cellId)")

        fig[:add_subplot](4,2,(i-1)*2+2)
        PyPlot.plt[:hist](synapseDistribution[i][2]./1000; range=postrange, bins=bins)
        PyPlot.plt[:xlabel]("distrance from postsynapse to soma ($cellId)")
        
    end 
    PyPlot.savefig(saveName)
end

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

function plot_connectivity_matrix(orderedConnMatrix::SparseMatrixCSC{T}; 
                                    title="connectivity matrix",
                                    max_size::Real=10) where T
    X,Y,V = SparseArrays.findnz(orderedConnMatrix)
    # transform to float and normalize to max_size
    V = Float32.(V)
    V ./= maximum(V) / Float32(max_size)

    trace = PlotlyJS.scatter(; x=X, y=Y, 
        mode="markers", 
        marker=attr(size=V, line_width=0, opacity=0.4))
    data = [trace,]

    layout = Layout(; title=title,
        width=800,
        height=800,
        xaxis = attr(title="presynapse", range=(1, size(orderedConnMatrix,1))),
        yaxis= attr(title="postsynapse", range=(1, size(orderedConnMatrix,2)))
    )
    return PlotlyJS.plot(data, layout)
end

function plot(neuron::Neuron; nodeStep::Integer=10, semantic::Bool=true, showSynapse::Bool=true)
    traces = PlotlyJS.GenericTrace[]
    # plot soma
    root = Neurons.get_root_node(neuron)
    root_trace = PlotlyJS.scatter3d(;x=[root[1]], y=[root[2]], z=[root[3]], 
                                    mode="markers", marker_size=7)
    push!(traces, root_trace)

    for (i,segment) in enumerate(neuron)
        nodeList = Neurons.get_segment_node_list(neuron, i)
        x = map(n->n[1], nodeList[1:nodeStep:end]) |> Vector{Float32}
        y = map(n->n[2], nodeList[1:nodeStep:end]) |> Vector{Float32}
        z = map(n->n[3], nodeList[1:nodeStep:end]) |> Vector{Float32}
        if semantic
            class = Neurons.Segments.get_class(segment)
            if class == Neurons.Segments.AXON_CLASS
                color = "rgb(255,0,0)"
            elseif class == Neurons.Segments.DENDRITE_CLASS
                color = "rgb(0,0,255)"
            else
                color = "rgb(127,127,127)"
            end
            segmentTrace = PlotlyJS.scatter3d(;x=x,y=y,z=z, mode="lines", 
                                                line_color=color)
        else
            segmentTrace = PlotlyJS.scatter3d(;x=x,y=y,z=z, mode="lines")
        end
        push!(traces, segmentTrace)

        if showSynapse
            preSynapseList = Neurons.Segments.concatenate_pre_synapses(segment)
            postSynapseList = Neurons.Segments.concatenate_post_synapses(segment)
            preTrace = PlotlyJS.scatter3d(; 
                                x=map(c->Synapses.get_psd_coordinate(c)[1], preSynapseList),
                                y=map(c->Synapses.get_psd_coordinate(c)[2], preSynapseList),
                                z=map(c->Synapses.get_psd_coordinate(c)[3], preSynapseList),
                                mode="markers", marker_color="rgb(179,66,244)", marker_size=2)
            postTrace = PlotlyJS.scatter3d(; 
                                x=map(c->Synapses.get_psd_coordinate(c)[1], postSynapseList),
                                y=map(c->Synapses.get_psd_coordinate(c)[2], postSynapseList),
                                z=map(c->Synapses.get_psd_coordinate(c)[3], postSynapseList),
                                mode="markers", marker_color="rgb(66,223,244)", marker_size=2)
            push!(traces, preTrace)
            push!(traces, postTrace)
        end
    end
    layout = PlotlyJS.Layout(; showlegend=false)
    PlotlyJS.plot(traces, layout)
end 

function plot_v2(neuron::Neuron; nodeStep::Integer=10)
    #using PyPlot
    #PyPlot.pygui(true)

    # plot soma
    root = Neurons.get_root_node(neuron)
    PyPlot.scatter3D([root[1]], [root[2]], [root[3]])
        
    for segment in neuron
        nodeList = Neurons.Segments.get_node_list(segment)
        x = map(n->n[1], nodeList[1:nodeStep:end])
        y = map(n->n[2], nodeList[1:nodeStep:end])
        z = map(n->n[3], nodeList[1:nodeStep:end])
        PyPlot.plot(x,y,z)
    end
end 

function plot_v1(neuron::Neuron; nodeStep::Integer=10)
    segmentList = neuron.segmentList
    plotly()
    for branch in segmentList
        nodeList = Neurons.Segments.get_node_list(branch)
        x = map(n->n[1], nodeList[1:nodeStep:end])
        y = map(n->n[2], nodeList[1:nodeStep:end])
        z = map(n->n[3], nodeList[1:nodeStep:end])
        plot!(x,y,z) #color=rand(Colors.RGB))
    end 

    root = Neurons.get_root_node(neuron)
    plot!([root[1]], [root[2]], [root[3]], m=(2, :circle), leg=false)
end 

"""
    plot_maximum_intensity_projection(vol::Array; colorMap="gray")
color map could be "jet" if you like colorful display. 
"""
function plot_maximum_intensity_projection(vol::Array{T,3}; colorMap="gray") where T
    fig = PyPlot.figure()
    fig[:add_subplot](2,2,1)
    xy = maximum(vol, 3)[:,:,1]
    PyPlot.imshow(xy, colorMap)
    # colorbar()
    fig[:add_subplot](2,2,2)
    xz = maximum(vol, 2)[:,1,:]
    PyPlot.imshow(xz, colorMap)
    fig[:add_subplot](2,2,3)
    yz = maximum(vol, 1)[1,:,:] |> rotl90
    PyPlot.imshow(yz, colorMap)
    PyPlot.colorbar()
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
    for i in 1: size(hc.merges, 1)
        x1 = pos[hc.merges[i,1]][1]
        x2 = pos[hc.merges[i,2]][1]
        append!(xs, [x1,x1,x2,x2])

        y1 = pos[hc.merges[i,1]][2]
        y2 = pos[hc.merges[i,2]][2]
        useheight ? h = hc.heights[i] : h = 1
        newy = maximum([y1,y2]) + h
        append!(ys, [y1,newy,newy,y2])
    end
    return (reshape(xs, 4, size(hc.merges, 1)), reshape(ys, 4, size(hc.merges, 1)))
end

function treepositions(hc::Hclust, useheight::Bool)
    order = indexmap(hc.order)
    positions = Dict{}()
    for (k,v) in order
        positions[-k] = (v, 0)
    end
    for i in 1:size(hc.merges,1)
        xpos = mean([positions[hc.merges[i,1]][1], positions[hc.merges[i,2]][1]])
        if hc.merges[i,1] < 0 && hc.merges[i,2] < 0
            useheight ? ypos = hc.heights[i] : ypos = 1
        else
            useheight ? h = hc.heights[i] : h = 1
            ypos = maximum([positions[hc.merges[i,1]][2], positions[hc.merges[i,2]][2]]) + h
        end

        positions[i] = (xpos, ypos)
    end
    return positions
end 

function plot(clust::Hclust)
	Plots.gr()
	#Plots.plotlyjs()
	Plots.plot(hclustplot(clust, true), seriestype=:path, color=:black,
    grid=false, legend=false) #,  xticks=classificationIdList[clust.order])
end 

end # module

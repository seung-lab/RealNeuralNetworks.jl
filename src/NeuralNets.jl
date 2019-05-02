module NeuralNets
using SparseArrays
using DataFrames
using LightGraphs
using MetaGraphs
#using Query
using ..Neurons 
using ..Neurons.Segments.Synapses
using ..SWCs

using Query

export NeuralNet 
const NeuralNet = MetaDiGraph 

"""
    NeuralNet(neuronId2neuron::Dict{Int, Neuron}; neuronIdList::Vector{Int})
construct neuron graph. 
Parameters:
    neuronId2neuron: the dictionary mapping ID to neuron 
    neuronIdList: the vertex id was ordered according to this list  
Return:
    net: the MetaDiGraph representing neural network.
"""
function NeuralNet(neuronId2neuron::Dict{Int, Neuron{T}}; 
                   neuronIdList::Vector{Int} = collect(keys(neuronId2neuron))) where T
    neuronIdSet = Set{Int}(neuronIdList)
    net = NeuralNet(length(neuronIdList))
    
    # name the vertices 
    neuronId2vertexId = Dict{Int,Int}()
    for (vertexId, neuronId) in enumerate(neuronIdList)
        neuronId2vertexId[neuronId] = vertexId
        neuron = neuronId2neuron[ neuronId ]
        set_props!(net, vertexId, Dict(:id=>neuronId, :skeleton=>neuron))
    end 

    # add edges 
    selfConnectNeuronIdSet = Set{Int}()
    for neuronId in neuronIdList 
        neuron = neuronId2neuron[ neuronId ]
        vertexId = neuronId2vertexId[ neuronId ]
        preSynapseList = Neurons.get_all_pre_synapse_list(neuron)
        for preSynapse in preSynapseList 
            postNeuronId = Synapses.get_post_synaptic_segmentation_id(preSynapse)
            if postNeuronId == neuronId 
                @warn("skipping self connection edge of neuron: ", neuronId)
                push!(selfConnectNeuronIdSet, neuronId)
            elseif postNeuronId in neuronIdSet 
                # this is an effective synaptic edge 
                postNeuronVertexId = neuronId2vertexId[ postNeuronId ]
                edge = (vertexId, postNeuronVertexId)
                add_edge!(net, edge)
                # get properties 
                if has_prop(net, edge[1], edge[2], :synapses)
                    prop = get_prop(net, edge[1], edge[2], :synapses)
                    totalPSDSize = get_prop(net, edge[1], edge[2], :total_psd_size)
                else 
                    prop = Vector{Synapse{T}}()
                    totalPSDSize = 0
                end 
                push!(prop, preSynapse)
                set_prop!(net, edge[1], edge[2], :synapses, prop)
                set_prop!(net, edge[1], edge[2], :num_synapses, length(prop))
                totalPSDSize += Synapses.get_psd_size(preSynapse)
                set_prop!(net, edge[1], edge[2], :total_psd_size, totalPSDSize)
            end
        end 
    end 
    @warn("skipped self connection edges from $(length(selfConnectNeuronIdSet)) neurons.")
    net 
end 

function NeuralNet( syn::DataFrame; neuronDict::Dict{Int, Neuron}=Dict{Int,Neuron}(), 
                                    cellIdList=get_cell_id_list(syn) )
    N = length(cellIdList)
    net = NeuralNet(N)
    # name the vertices
    for (i,cellId) in enumerate(cellIdList)
        if haskey(neuronDict, cellId)
            set_props!(net, i, Dict(:id=>cellId, :skeleton=>neuronDict[cellId]))
        else 
            @warn("no skeleton attatched: ", cellId)
        end 
    end 

    # filter out the valid set
    cellIdSet = Set(cellIdList)
	syn =  @from i in syn begin
        @where 	i.presyn_segid in cellIdSet && 
                i.postsyn_segid in cellIdSet
        @select i 
        @collect DataFrame
    end
    # construct a look up table
    cellId2vertexIdMap = Dict{Int, Int}()
    for (i,cellId) in enumerate(cellIdList)
        cellId2vertexIdMap[cellId] = i
    end 
    # add edges
    numSelfConnection = 0
    for row in eachrow(syn)
        preId  = round(Int, row[:presyn_segid])
        postId = round(Int, row[:postsyn_segid])
        if preId == postId 
            numSelfConnection += 1
            @warn("self connection synapse on neuron: $(preId), skipping the connection")
            continue
        end 
        edge = ( cellId2vertexIdMap[preId], cellId2vertexIdMap[postId] )
        add_edge!(net, edge)
        # get properties
        if has_prop(net, edge[1], edge[2], :synapses)
            prop = get_prop(net, edge[1], edge[2], :synapses)
        else 
            prop = Vector{Synapse}()
        end 
        synapse = Synapse(row)
        push!(prop, synapse )
        set_prop!(net, edge[1], edge[2], :synapses, prop)
    end
    println("number of self connection synapses: $(numSelfConnection)")
    net 
end

function get_cell_id_list(syn::DataFrame)
    cellIdSet = Set(syn[:presyn_segid]) ∪ Set(syn[:postsyn_segid]) 
    return [map(x->round(Int, x), cellIdSet)...]
end 

function get_cell_id_list(self::NeuralNet)
    N = nv(self)
    cellIdList = Vector{Int}(N)
    for v in vertices(self)
        push!(cellIdList, get_prop(self, v, :id))
    end 
    cellIdList 
end 

"""
    get_synapse_number_connectivity_matrix(self::NeuralNet)
Return: 
    sparse matrix representing the connectivity matrix. 
    The rows are presynase and the columns are post synapses.
    For example, [i,j] means the axon of neuron i innervate the dendrite of neuron j
"""
function get_synapse_number_connectivity_matrix(self::NeuralNet)
    N = nv(self)
    ret = spzeros(Int, N,N)
    for edge in edges(self)
        synapseList = get_prop(self, edge, :synapses)
        ret[src(edge), dst(edge)] = length(synapseList)
    end 
    ret 
end

"""
    get_psd_size_connectivity_matrix(self::NeuralNet)
Return: 
    sparse matrix representing the connectivity matrix.
    The rows are presynapse and the columns are post synapses.
    For example [i,j] means the axon of neuron i innvervate the dendrite of neuron j.
    The value represent the total volume of psd
"""
function get_psd_size_connectivity_matrix(self::NeuralNet)
    N = nv(self)
    ret = spzeros(Int, N,N)
    for edge in edges(self)
        synapses = get_prop(self, edge, :synapses)
        totalSynSize = 0
        for synapse in synapses
            # note that the size was measured using voxel number
            totalSynSize += Synapses.get_psd_size(synapse)
        end 
        ret[src(edge), dst(edge)] = totalSynSize 
    end 
    ret 
end 

"""
    labels2orders(label::Vector)

transform clustering label to orders 
    
Parameters:
    labels: the assignment of each neuron. The label must start from 1 to N, where N is the maximum number of groups. 
"""
function labels2orders(labels::Vector)
    orders = zero(labels)
    i = 0
    @inbounds for l in unique(labels)
        idxes = findall(x->x==l, labels)
        for idx in idxes
            i += 1
            orders[idx] = i
        end
    end 
    return orders
end 

"""
reorder_connectivity_matrix(connMatrix::Union{SparseMatrixCSC{T}, Matrix{T}}, orders::Vector) 

reorder the rows and columns of the connectivity matrix for visualization.

Parameters:
    connMatrix: the connectivity matrix
    orders: the rearranged order of samples. 

Return:
    a reordered connectivity matrix which is also a SparseMatrixCSC 
"""
function reorder_matrix(connMatrix::Union{SparseMatrixCSC{T},Matrix{T}}, orders::Vector) where T 

    # the reordered connectivity Matrix
    if isa(connMatrix, SparseMatrixCSC)
        ret = spzeros(T, size(connMatrix)...,)
    else 
        ret = zeros(T, size(connMatrix))
    end
    
    @inbounds for (i, o1) in enumerate(orders)
        for (j, o2) in enumerate(orders)
            if connMatrix[o1,o2] !== zero(T)
                ret[i,j] = connMatrix[o1,o2]
            end 
        end
    end 
    return ret
end 

function remove_singletons(self::NeuralNet)
    @error("unimplemented")
end 


function get_nblast_similarity_matrix(self::NeuralNet)
    error("not implemented.")
end 



end # end of module

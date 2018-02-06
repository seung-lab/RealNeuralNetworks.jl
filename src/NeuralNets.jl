module NeuralNets
using DataFrames
using LightGraphs
using MetaGraphs
#using Query
using ..Neurons 

export NeuralNet 
const NeuralNet = MetaDiGraph 

function get_cell_id_list(syn::DataFrame)
    cellIdSet = Set(syn[:presyn_seg]) ∪ Set(syn[:postsyn_seg]) 
    return [map(x->round(Int, x), cellIdSet)...]
end 

function NeuralNet( neuronDict::Dict{Int, Neuron}, syn::DataFrame; 
                                        cellIdList=get_cell_id_list(syn) )
    const N = length(cellIdList)
    net = NeuralNet(N)
    # name the vertices
    for (i,cellId) in enumerate(cellIdList) 
        set_props!(net, i, Dict(:id=>cellId, :skeleton=>neuronDict[cellId]))
    end 

#    # filter out the valid set 
#	syn =  @from i in syn begin
#        @where 	i.presyn_seg  <=maximum(cellIdList) && 
#                i.presyn_seg  >=minimum(cellIdList) && 
#				i.postsyn_seg <=maximum(cellIdList) && 
#                i.postsyn_seg >=minimum(cellIdList) &&
#           	    round(Int,i.presyn_seg.value) in cellIdSet && 
#                round(Int,i.postsyn_seg.value) in cellIdSet
#        @select i 
#        @collect DataFrame
#    end
    # construct a look up table
    cellId2vertexIdMap = Dict{Int, Int}()
    for (i,cellId) in enumerate(cellIdList)
        cellId2vertexIdMap[cellId] = i
    end 
    # add edges
    numSelfConnection = 0
    for row in eachrow(syn)
        preId  = round(Int, row[:presyn_seg])
        postId = round(Int, row[:postsyn_seg])
        if preId == postId 
            numSelfConnection += 1
            println("self connection synapse on neuron: $(preId), skipping the connection")
            continue
        end 
        edge = ( cellId2vertexIdMap[preId], cellId2vertexIdMap[postId] )
        add_edge!(net, edge)
        # get properties
        if has_prop(net, edge[1], edge[2], :synapses)
            prop = get_prop(net, edge[1], edge[2], :synapses)
        else 
            prop = Vector{NTuple{23,Int}}()
        end 
        push!(prop, (map(x->round(Int,x[2]), row)...) )
        set_prop!(net, edge[1], edge[2], :synapses, prop)
    end
    println("number of self connection synapses: $(numSelfConnection)")
    net 
end

function get_cell_id_list(self::NeuralNet)
    N = nv(self)
    cellIdList = Vector{Int}(N)
    for v in vertices(self)
        push!(cellIdList, get_prop(self, v, :id))
    end 
    cellIdList 
end 

function get_synapse_number_connectivity_matrix(self::NeuralNet)
    N = nv(self)
    ret = spzeros(Int, N,N)
    for edge in edges(self)
        ret[src(edge), dst(edge)] = length(get_prop(self, edge, :synapses))
    end 
    ret 
end

function get_synapse_size_connectivity_matrix(self::NeuralNet)
    N = nv(self)
    ret = spzeros(Int, N,N)
    for edge in edges(self)
        synapses = get_prop(self, edge, :synapses)
        totalSynSize = 0
        for synapse in synapses 
            totalSynSize += synapse[23]
        end 
        ret[src(edge), dst(edge)] = totalSynSize 
    end 
    ret 
end 

function synapse_num_weighted!(self::NeuralNet)
    error("unimplemented")
end 

end # end of module

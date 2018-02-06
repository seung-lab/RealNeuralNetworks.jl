module SynapseTables
using DataFrames
using Query
using ..Utils.BoundingBoxes

export SynapseTable 
const SynapseTable = DataFrame  

function preprocessing!(self::SynapseTable)
	DataFrames.dropmissing!(self)
	for (key, value) in DataFrames.eachcol(self)
		if key!=:presyn_wt && key!=:postsyn_wt
			self[key] = round.(Int, value)
		end
	end
    self[:BBOX_bx] .+= 1
	self[:BBOX_by] .+= 1
	self[:BBOX_bz] .+= 1
	self[:BBOX_ex] .+= 1
	self[:BBOX_ey] .+= 1
	self[:BBOX_ez] .+= 1
	self[:COM_x]   .+= 1
	self[:COM_y]   .+= 1
	self[:COM_z]   .+= 1
	self[:presyn_x].+= 1
	self[:presyn_y].+= 1
	self[:presyn_z].+= 1
	self[:postsyn_x].+= 1
	self[:postsyn_y].+= 1
	self[:postsyn_z].+= 1
	nothing
end 

function get_synaptic_boutons_of_a_neuron(self::SynapseTable, neuronId::Int)
    ret =  @from i in self begin                                         
        @where i.presyn_seg == neuronId 
        @select i                                                      
        @collect DataFrame                                             
    end                                                                
    ret
end 

function get_postsynaptic_density_of_a_neuron(self::SynapseTable, neuronId::Int)
    ret =  @from i in self begin                                         
        @where i.postsyn_seg == neuronId
        @select i                                                      
        @collect DataFrame                                             
    end                                                                
    ret
end 

function get_synapses_of_a_neuron(self::SynapseTable, neuronId::Int)
	ret =  @from i in self begin                                         
        @where i.postsyn_seg == neuronId || i.presyn_seg == neuronId
        @select i                                                      
        @collect DataFrame                                             
    end                                                                
    ret
end 

function BoundingBox(self::SynapseTable)
    minCorner = ( minimum(self[:COM_x]), minimum(self[:COM_y]), minimum(self[:COM_z]))
    maxCorner = ( maximum(self[:COM_x]), maximum(self[:COM_y]), maximum(self[:COM_z]))
    BoundingBoxes.BoundingBox(minCorner, maxCorner)
end 

end # module

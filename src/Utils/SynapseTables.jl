module SynapseTables
using DataFrames
using Query
using ..Utils.BoundingBoxes
using OffsetArrays

export SynapseTable 
const SynapseTable = DataFrame  

function preprocessing!(self::SynapseTable, voxelSize::NTuple{3,Int})
	DataFrames.dropmissing!(self)
	for (key, value) in DataFrames.eachcol(self)
		if key!=:presyn_wt && key!=:postsyn_wt
			self[key] = round.(Int, value)
		end
    end
    
    DataFrames.rename!(self, [  :centroid_x => :psd_x, 
                                :centroid_y => :psd_y,
                                :centroid_z => :psd_z])

    # transform to physical coordinate
    self[:BBOX_bx]      .*= voxelSize[1]
    self[:BBOX_by]      .*= voxelSize[2]
    self[:BBOX_bz]      .*= voxelSize[3]
    self[:BBOX_ex]      .*= voxelSize[1]
    self[:BBOX_ey]      .*= voxelSize[2]
    self[:BBOX_ez]      .*= voxelSize[3]
    self[:psd_x]        .*= voxelSize[1]
    self[:psd_y]        .*= voxelSize[2]
    self[:psd_z]        .*= voxelSize[3]
    self[:presyn_x]     .*= voxelSize[1]
    self[:presyn_y]     .*= voxelSize[2]
    self[:presyn_z]     .*= voxelSize[3]
    self[:postsyn_x]    .*= voxelSize[1]
    self[:postsyn_y]    .*= voxelSize[2]
    self[:postsyn_z]    .*= voxelSize[3]
	nothing
end

function postprocessing!(self::SynapseTable, voxelSize::NTuple{3,Int})
    # transform to voxel coordinate
    self[:BBOX_bx]      ./= voxelSize[1]
    self[:BBOX_by]      ./= voxelSize[2]
    self[:BBOX_bz]      ./= voxelSize[3]
    self[:BBOX_ex]      ./= voxelSize[1]
    self[:BBOX_ey]      ./= voxelSize[2]
    self[:BBOX_ez]      ./= voxelSize[3]
    self[:psd_x]        ./= voxelSize[1]
    self[:psd_y]        ./= voxelSize[2]
    self[:psd_z]        ./= voxelSize[3]
    self[:presyn_x]     ./= voxelSize[1]
    self[:presyn_y]     ./= voxelSize[2]
    self[:presyn_z]     ./= voxelSize[3]
    self[:postsyn_x]    ./= voxelSize[1]
    self[:postsyn_y]    ./= voxelSize[2]
    self[:postsyn_z]    ./= voxelSize[3]
end 

function get_coordinate_array(self::SynapseTable, prefix::String)
	coordinateNames = get_coordinate_names(prefix)
	hcat(	self[coordinateNames[1]], 
			self[coordinateNames[2]], 
			self[coordinateNames[3]])
end 

function get_coordinate(self::DataFrameRow, prefix::String)
	coordinateNames = get_coordinate_names(prefix)
    map(x->self[x], coordinateNames)
end 


function get_pre_synapses_of_a_neuron(self::SynapseTable, neuronId::Int)
    ret =  @from i in self begin                                         
        @where i.presyn_segid == neuronId 
        @select i                                                      
        @collect DataFrame                                             
    end                                                                
    ret
end 

function get_post_synapses_of_a_neuron(self::SynapseTable, neuronId::Int)
    ret =  @from i in self begin                                         
        @where i.postsyn_segid == neuronId 
        @select i                                                      
        @collect DataFrame                                             
    end                                                                
    ret
end 


function get_synapses_of_a_neuron(self::SynapseTable, neuronId::Int)
	ret =  @from i in self begin                                         
        @where i.postsyn_segid == neuronId || i.presyn_segid == neuronId
        @select i                                                      
        @collect DataFrame                                             
    end                                                                
    ret
end 

@inline function get_coordinate_names( prefix::String )
    map(x->Symbol(prefix*x), ("x", "y", "z"))
end 

function initialize_mask(self::SynapseTable, voxelSize::NTuple{3,Int};
                coordinatePrefixList = ["presyn_", "postsyn_"], T::DataType = Bool)
	@assert !isempty(self)
    range = mapreduce(x->BoundingBox(self,x), union, coordinatePrefixList) |> BoundingBoxes.get_unit_range
    range = map((r,s)->fld(r.start, s):cld(r.stop, s), range, voxelSize)
    mask = OffsetArray{T}(range...)
    fill!(mask, zero(T))
    mask
end

function get_mask(self::SynapseTable; 
                voxelSize::NTuple{3,Int} = (1000, 1000, 1000) #=nm=#,
                coordinatePrefixList::Vector{String} = ["presyn_", "postsyn_"],
                mask::OffsetArray{T,3,Array{T,3}} = initialize_mask(self, voxelSize; 
                    coordinatePrefixList = coordinatePrefixList, T=Float32)) where T
    @assert !isempty(self)
	for coordinatePrefix in coordinatePrefixList 
		coordinateNames = SynapseTables.get_coordinate_names( coordinatePrefix )
		coordinateList = map((c,s)->div.(self[c], s), coordinateNames, voxelSize)
		for i in 1:length(coordinateList[1])
			mask[ 	coordinateList[1][i], 
					coordinateList[2][i], 
					coordinateList[3][i] ] = one(eltype(mask))
		end
	end
    mask
end

@inline function BoundingBox(self::SynapseTable, prefix::String)
    coordinateNames = get_coordinate_names( prefix )
    BoundingBox(self, coordinateNames)
end 

@inline function BoundingBox(self::SynapseTable, coordinateNames::NTuple{3,Symbol})
    if isempty(self)
        warn("empty synapse table!")
        return nothing
    else 
        minCorner = map(x->minimum(self[x]), coordinateNames)
        maxCorner = map(x->maximum(self[x]), coordinateNames)
        return BoundingBoxes.BoundingBox(minCorner, maxCorner)
    end 
end 



end # module

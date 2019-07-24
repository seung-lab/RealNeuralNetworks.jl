module SynapseTables
using DataFrames
using Query
using ..Utils.BoundingBoxes
using OffsetArrays
using ProgressMeter 

export SynapseTable 
const SynapseTable = DataFrame  

function preprocess(self::SynapseTable, voxelSize::Tuple)
    @assert length(voxelSize) == 3
    voxelSize = map(Float32, voxelSize)
    
    println("drop missing values...")
	@time DataFrames.dropmissing!(self)

    # this is not working due to a DataFrames bug 
    # https://github.com/JuliaData/DataFrames.jl/issues/1862
    # we need to add true to make it work as old ways...
	@showprogress 1 "transform datatype to Float32..." for (key, value) in DataFrames.eachcol(self, true)
		#if key!=:presyn_wt && key!=:postsyn_wt
		#	self[key] = round.(Int, value)
        #end
        self[!, key] = Vector{Float32}(value)
    end

    println("remove self connections...")
    nrow1 = DataFrames.nrow(self)
    self = @from i in self begin 
        @where i.presyn_segid != i.postsyn_segid 
        @select i 
        @collect DataFrame 
    end  
    nrow2 = DataFrames.nrow(self)
    println("removed $(nrow1-nrow2) self connections in $nrow1 synaptic connections.")
   
    try 
        DataFrames.rename!(self, [  :centroid_x => :psd_x, 
                                    :centroid_y => :psd_y,
                                    :centroid_z => :psd_z])
    catch err 
        @info("did not find centroid column in the synapse table file.")
    end

    println("transform to physical coordinate...")
    self[!, :BBOX_bx]      = (self[!, :BBOX_bx]   .- one(Float32)) .* voxelSize[1]
    self[!, :BBOX_by]      = (self[!, :BBOX_by]   .- one(Float32)) .* voxelSize[2]
    self[!, :BBOX_bz]      = (self[!, :BBOX_bz]   .- one(Float32)) .* voxelSize[3]
    
    self[!, :BBOX_ex]      = (self[!, :BBOX_ex]   .- one(Float32)) .* voxelSize[1]
    self[!, :BBOX_ey]      = (self[!, :BBOX_ey]   .- one(Float32)) .* voxelSize[2]
    self[!, :BBOX_ez]      = (self[!, :BBOX_ez]   .- one(Float32)) .* voxelSize[3]

    self[!, :psd_x]        = (self[!, :psd_x]     .-one(Float32)) .* voxelSize[1]
    self[!, :psd_y]        = (self[!, :psd_y]     .-one(Float32)) .* voxelSize[2]
    self[!, :psd_z]        = (self[!, :psd_z]     .-one(Float32)) .* voxelSize[3]
    
    self[!, :presyn_x]     = (self[!, :presyn_x]  .- one(Float32)) .* voxelSize[1]
    self[!, :presyn_y]     = (self[!, :presyn_y]  .- one(Float32)) .* voxelSize[2]
    self[!, :presyn_z]     = (self[!, :presyn_z]  .- one(Float32)) .* voxelSize[3]
    
    self[!, :postsyn_x]    = (self[!, :postsyn_x] .- one(Float32)) .* voxelSize[1]
    self[!, :postsyn_y]    = (self[!, :postsyn_y] .- one(Float32)) .* voxelSize[2]
    self[!, :postsyn_z]    = (self[!, :postsyn_z] .- one(Float32)) .* voxelSize[3]
	return self
end

function postprocess(self::SynapseTable, voxelSize::Union{Vector, Tuple})
    @assert length(voxelSize) == 3
    voxelSize = map(Float32, voxelSize)

    println("transform to voxel coordinate...")
    self[!, :BBOX_bx]      = (self[!, :BBOX_bx]./voxelSize[1]) .+ one(Float32)
    self[!, :BBOX_by]      = (self[!, :BBOX_by]./voxelSize[2]) .+ one(Float32)
    self[!, :BBOX_bz]      = (self[!, :BBOX_bz]./voxelSize[3]) .+ one(Float32)
    
    self[!, :BBOX_ex]      = (self[!, :BBOX_ex]./voxelSize[1]) .+ one(Float32)
    self[!, :BBOX_ey]      = (self[!, :BBOX_ey]./voxelSize[2]) .+ one(Float32)
    self[!, :BBOX_ez]      = (self[!, :BBOX_ez]./voxelSize[3]) .+ one(Float32)

    
    self[!, :psd_x]        = (self[!, :psd_x]./voxelSize[1]) .+ one(Float32)
    self[!, :psd_y]        = (self[!, :psd_y]./voxelSize[2]) .+ one(Float32)
    self[!, :psd_z]        = (self[!, :psd_z]./voxelSize[3]) .+ one(Float32)
    
    self[!, :presyn_x]     = (self[!, :presyn_x]./voxelSize[1]) .+ one(Float32)
    self[!, :presyn_y]     = (self[!, :presyn_y]./voxelSize[2]) .+ one(Float32)
    self[!, :presyn_z]     = (self[!, :presyn_z]./voxelSize[3]) .+ one(Float32)
 
    self[!, :postsyn_x]     = (self[!, :postsyn_x]./voxelSize[1]) .+ one(Float32)
    self[!, :postsyn_y]     = (self[!, :postsyn_y]./voxelSize[2]) .+ one(Float32)
    self[!, :postsyn_z]     = (self[!, :postsyn_z]./voxelSize[3]) .+ one(Float32)
    
    return self
end 

function get_coordinate_array(self::SynapseTable, prefix::String)
	coordinateNames = get_coordinate_names(prefix)
	hcat(	self[!, coordinateNames[1]], 
			self[!, coordinateNames[2]], 
			self[!, coordinateNames[3]])
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
    mask = OffsetArray{T}(undef, range...)
    fill!(mask, zero(T))
    mask
end

function get_mask(self::SynapseTable; 
                  voxelSize::Union{Vector, Tuple} = (1000, 1000, 1000) #=nm=#,
                coordinatePrefixList::Vector{String} = ["presyn_", "postsyn_"],
                mask::OffsetArray{T,3,Array{T,3}} = initialize_mask(self, voxelSize; 
                    coordinatePrefixList = coordinatePrefixList, T=Float32)) where T
    @assert !isempty(self)
    @assert length(voxelSize) == 3
    voxelSize = map(Float32, voxelSize)
	for coordinatePrefix in coordinatePrefixList 
		coordinateNames = SynapseTables.get_coordinate_names( coordinatePrefix )
		coordinateList = map((c,s)->self[!, c]/s, coordinateNames, voxelSize)
        coordinateList = map(x->round.(Int, x), coordinateList)
		for i in 1:length(coordinateList[1])
			mask[   coordinateList[1][i], 
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
        @warn("empty synapse table!")
        return nothing
    else 
        minCorner = map(x->minimum(self[!, x]), coordinateNames)
        maxCorner = map(x->maximum(self[!, x]), coordinateNames)
        return BoundingBoxes.BoundingBox(minCorner, maxCorner)
    end 
end 



end # module

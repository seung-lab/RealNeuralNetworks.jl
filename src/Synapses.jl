module Synapses
using RealNeuralNetworks.Utils.SynapseTables 
using RealNeuralNetworks.Utils.BoundingBoxes
using DataFrames

export Synapse 

mutable struct Synapse{T}
    psdSegmentationId               ::Int 
    psdCoordinate                   ::NTuple{3,T} # nm
    psdBoundingBox                  ::BoundingBox 
    psdSize                         ::Int 
    preSynapticSegmentationId       ::Int 
    preSynapticCoordinate           ::NTuple{3,T} # nm
    preSynapticWeight               ::Float32 
    postSynapticSegmentationId      ::Int  
    postSynapticCoordinate          ::NTuple{3,T} # nm
    postSynapticWeight              ::Float32 
end 

function Synapse( row::DataFrameRow )
    psdBoundingBox = BoundingBox(   SynapseTables.get_coordinate(row, "BBOX_b"),
                                    SynapseTables.get_coordinate(row, "BBOX_e") )
    psdCoordinate = SynapseTables.get_coordinate(row, "psd_")
    preSynapticCoordinate = SynapseTables.get_coordinate(row, "presyn_") 
    postSynapticCoordinate = SynapseTables.get_coordinate(row, "postsyn_")

    Synapse{Float32}(row[:psd_segid], psdCoordinate, psdBoundingBox, row[:size], 
            row[:presyn_segid],  preSynapticCoordinate,  row[:presyn_wt],
            row[:postsyn_segid], postSynapticCoordinate, row[:postsyn_wt])
end

function Synapse( df::DataFrame )
    @assert DataFrames.nrow( df ) == 1
    Synapse( df[1,:] )
end 

############### Base functions ################
function Base.show(self::Synapse)
    println("psd segmentation id: ", psdSegmentationId)
    println("psd coordinate: ", psdCoordinate)
    println("psd bounding box: ", psdBoundingBox)
    println("psd size: ", psdSize)
    println("\npre synaptic segmentation id: ", preSynapticSegmentationId)
    println("pre synaptic coordinate: ", preSynapticCoordinate)
    println("pre synaptic weight: ", preSynapticWeight)
    println("\npost synaptic segmentation id: ", postSynapticSegmentationId)
    println("post synaptic coordinate: ", postSynapticCoordinate)
    println("post synaptic weight: ", postSynapticWeight)
    nothing
end 

@inline function Base.zero(self::Synapse) 
    nothing 
end

############### properties ####################
"""
Note that the psd size was measured using voxel numbers.
"""
@inline function get_psd_size(self::Synapse) self.psdSize end  
@inline function get_psd_segmentation_id(self::Synapse) self.psdSegmentationId end 
@inline function get_psd_coordinate( self::Synapse ) self.psdCoordinate end 
@inline function get_psd_bounding_box( self::Synapse ) self.psdBoundingBox end 

@inline function get_pre_synaptic_segmentation_id( self::Synapse ) self.preSynapticSegmentationId end 
@inline function get_pre_synaptic_coordinate( self::Synapse ) self.preSynapticCoordinate end
@inline function get_pre_synaptic_weight( self::Synapse ) self.preSynapticWeight end 

@inline function get_post_synaptic_segmentation_id( self::Synapse ) 
    self.postSynapticSegmentationId 
end 
@inline function get_post_synaptic_coordinate( self::Synapse ) 
    self.postSynapticCoordinate 
end
@inline function get_post_synaptic_weight( self::Synapse ) self.postSynapticWeight end 

#@inline function isbutton(self::Synapse, neuronId::Int)
#    neuronId == get_presynaptic_segmentation_id(self)
#end 

end # module

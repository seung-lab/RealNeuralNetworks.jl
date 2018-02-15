module Synapses
using RealNeuralNetworks.Utils.SynapseTables 
using RealNeuralNetworks.Utils.BoundingBoxes
using DataFrames

export Synapse 

struct Synapse
    psdSegmentationId               ::Int 
    psdCoordinate                   ::NTuple{3,Float32} # nm
    psdBoundingBox                  ::BoundingBox 
    preSynapticSegmentationId       ::Int 
    preSynapticCoordinate           ::NTuple{3,Float32} # nm
    preSynapticWeight               ::Float32 
    postSynapticSegmentationId      ::Int 
    postSynapticCoordinate          ::NTuple{3,Float32} # nm
    postSynapticWeight              ::Float32 
end 

function Synapse( row::DataFrameRow )
    psdBoundingBox = BoundingBox(   SynapseTables.get_coordinate(row, "BBOX_b"),
                                    SynapseTables.get_coordinate(row, "BBOX_e") )
    psdCoordinate = SynapseTables.get_coordinate(row, "psd_")
    preSynapticCoordinate = SynapseTables.get_coordinate(row, "presyn_") 
    postSynapticCoordinate = SynapseTables.get_coordinate(row, "postsyn_")

    Synapse(row[:psd_segid], psdCoordinate, psdBoundingBox, 
            row[:presyn_segid],  preSynapticCoordinate,  row[:presyn_wt],
            row[:postsyn_segid], postSynapticCoordinate, row[:postsyn_wt])
end

function Synapse( df::DataFrame )
    @assert DataFrames.nrow( df ) == 1
    Synapse( df[1,:] )
end 

############### properties ####################
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


@inline function isbutton(self::Synapse, neuronId::Int)
    neuronId == get_presynaptic_segmentation_id(self)
end 

end # module

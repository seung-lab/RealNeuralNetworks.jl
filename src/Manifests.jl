module Manifests
using GSDicts, BigArrays
export Manifest

immutable Manifest
    # the bigarray for cutout
    ba          ::AbstractBigArray
    # the id of object
    objId       ::Integer
    # the unit range list for cutout
    rangeList   ::Tuple{UnitRange}
end 

"""
Parameters:
    dir 
"""
function Manifest( manifestDirPath::AbstractString, manifestKey::AbstractString, bigArrayPath::AbstractString )
    ba = BigArray( GSDict( bigArrayPath ) )
    h = GSDicts( manifestDirPath; valueType=Dict{Symbol, Any})
    Manifest( h[manifestKey], ba )
end
"""
example: {"fragments": ["770048087:0:2968-3480_1776-2288_16912-17424"]}
"""
function Manifest( h::Dict{Symbol, Any}, ba::AbstractBigArray )
    Manifest( h[:fragments], ba )
end 
"""
example: ["770048087:0:2968-3480_1776-2288_16912-17424"]
"""
function Manifest{D,T,N,C}( ranges::Vector{String}, ba{D,T,N,C} )
    objId = parse( split(ranges[1], ":")[1] )
    objId = convert(T, objId)
    ranges = map(x-> split(x,":")[end], ranges)
    rangeList = map( BigArrays.Indexes.string2unit_range, ranges )
    Manifest( ba, objId, rangeList )
end

function Base.start(self::Manifest)
    1
end 

"""
get the point cloud and dbf
"""
function Base.next(self::Manifest, i )
    range = Manifest.rangeList[i]
    seg = self.ba[range ...]
     
    return (point_cloud, dbf), i+1
end 
function Base.done(self::Manifest, i)
    i > length( self.rangeList )
end 

end # module

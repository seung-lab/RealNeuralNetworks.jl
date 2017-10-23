module PointObjs
export PointObj

type PointObj
    point_type  :: UInt8 
    x           :: Float32 
    y           :: Float32 
    z           :: Float32 
    radius      :: Float32 
    parent      :: Int32
end

function PointObj( p::Union{Tuple, Vector} )
    @assert length(p)==6
    PointObj( UInt8(p[1]), Float32(p[2]), Float32(p[3]), Float32(p[4]), Float32(p[5]), Int32(p[6]) )
end 

function get_parent( self::PointObj )
    self.parent 
end 

function Base.String(self::PointObj)
    "$(self.point_type) $(self.x) $(self.y) $(self.z) $(self.radius) $(self.parent)"
end 

function euclidean_distance(self::PointObj, other::PointObj) 
    norm( [self.x - other.x, self.y - other.y, self.z - other.z] )
end 

end # end of module

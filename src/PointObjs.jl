module PointObjs
export PointObj

mutable struct PointObj
    class       :: UInt8 
    x           :: Float32 
    y           :: Float32 
    z           :: Float32 
    radius      :: Float32 
    parent      :: Int32
end

function PointObj( p::Tuple )
    @assert length(p) == 6
    PointObj( UInt8(p[1]), Float32(p[2]), Float32(p[3]), Float32(p[4]), Float32(p[5]), Int32(p[6]) )
end

function PointObj( data::Vector{UInt8} )
    @assert length(data) == 21
    PointObj(   data[1], 
             reinterpret(Float32, data[2:5])[1],
             reinterpret(Float32, data[6:9])[1],
             reinterpret(Float32, data[10:13])[1],
             reinterpret(Float32, data[14:17])[1],
             reinterpret(Int32, data[18:21])[1])
end 

function get_parent( self::PointObj )
    self.parent 
end 

function serialize(self::PointObj)
    io = IOBuffer(read=false, write=true, maxsize=21)
    write(io, self.class)
    write(io, self.x)
    write(io, self.y)
    write(io, self.z)
    write(io, self.radius)
    write(io, self.parent)
    take!(io)
end 

function Base.String(self::PointObj)
    "$(self.class) $(self.x) $(self.y) $(self.z) $(self.radius) $(self.parent)"
end

function Base.isequal(self::PointObj, other::PointObj)
    self.class == other.class && 
    self.x == other.x &&
    self.y == other.y &&
    self.z == other.z &&
    self.parent == other.parent 
end 

function Base.:(==)(self::PointObj, other::PointObj) isequal(self, other) end 

function euclidean_distance(self::PointObj, other::PointObj) 
    norm( [self.x - other.x, self.y - other.y, self.z - other.z] )
end 

end # end of module

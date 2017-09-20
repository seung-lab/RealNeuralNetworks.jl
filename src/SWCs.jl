module SWCs
export SWC

type PointObj
    point_type  :: UInt8 
    x           :: Float32 
    y           :: Float32 
    z           :: Float32 
    radius      :: Float32 
    parent      :: Int32
end 

function to_string(self::PointObj)
    "$(self.point_type) $(self.x) $(self.y) $(self.z) $(self.radius) $(self.parent)"
end 

type SWC
    points      :: Vector{PointObj}
    attributes  :: Dict{Symbol, Any}
end

function SWC{T}(nodes::Array{T,2}, edges::Vector, 
             roots::Vector, radii::Vector, destinations::Vector; 
             attributes = Dict{Symbol,Any}())
    @assert size(nodes, 1) == length(radii)
    num = length(radii)
    points = Vector{PointObj}()
    sizehint!(points, num)

    for i in 1:num 
        point = PointObj(0, nodes[i,1], nodes[i,2], nodes[i,3], radii[i], -1)
        push!(points, point)
    end
    # assign parents according to edge 
    for e in edges 
        points[e[2]].parent = e[1]
    end  
    SWC(points, attributes)
end 

function get_points_num(self::SWC)  length(self.points) end 
function get_points(self::SWC) self.points end 

function save(self::SWC, file_name::AbstractString)
    points = get_points(self)
    f = open(file_name, "w")
    for i in 1:get_points_num(self)
        write(f, "$i $(to_string(points[i])) \n")
    end
    close(f)
end 

end # module of SWCs

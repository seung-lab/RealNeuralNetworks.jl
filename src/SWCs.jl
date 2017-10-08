module SWCs
using ..TEASAR.NodeNets

export SWC

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

function to_string(self::PointObj)
    "$(self.point_type) $(self.x) $(self.y) $(self.z) $(self.radius) $(self.parent)"
end 

type SWC
    # note that the first integer label was implicitly encoded as index
    points      :: Vector{PointObj}
end

function SWC()
    SWC( Vector{PointObj}() )
end 

function SWC(nodeNet::NodeNet)
    edges = NodeNets.get_edges(nodeNet)
    points = Vector{PointObj}()
    sizehint!(points, NodeNets.get_node_num(nodeNet))

    for node in NodeNets.get_nodes(nodeNet)
        point = PointObj(0, node[1], node[2], node[3], node[4], -1)
        push!(points, point)
    end
    # assign parents according to edge 
    for e in edges 
        points[e[2]].parent = e[1]
    end  
    SWC(points)
end 

################### get properties ###################
function Base.length(self::SWC)  length(self.points) end 
function get_points(self::SWC) self.points end 

################## IO ###############################
function save(self::SWC, file_name::AbstractString)
    points = get_points(self)
    f = open(file_name, "w")
    for i in 1:length(self)
        write(f, "$i $(to_string(points[i])) \n")
    end
    close(f)
end 

function load(file_name::AbstractString)
    swc = SWC()
    open(file_name) do f
        for line in eachline(f)
            try 
                numbers = map(parse, split(line))
                # construct a point object
                pointObj = PointObj( numbers[2:7] )
                push!(swc, pointObj)
            catch err 
                if !constains(line, "#")
                    println("comment in swc file: $line")
                else
                    warn("invalid line: $line")
                end 
            end 
        end 
    end 
    return swc
end 

#################### manipulate ######################
"""
add more point object 
"""
function Base.push!(self::SWC, p::PointObj )
    push!(self.points, p)
end 

"""
note that only stretch the coordinates here, not including the radius
since radius was already adjusted in the neighborhood weights 
"""
function stretch_coordinates!(self::SWC, expansion::Tuple)
    @assert length(expansion) == 3
    for i in 1:length( self.points )
        self.points[i].x *= expansion[1]
        self.points[i].y *= expansion[2]
        self.points[i].z *= expansion[3]
        self.points[i].radius *= (prod(expansion))^(1/3)
    end 
end 

"""
stretch the coordinate according to the mip level
normally, we only build mip level at XY plane, not Z
"""
function stretch_coordinates!(self::SWC, mip::Integer)
    stretch_coordinates!(self, (2^mip, 2^mip, 1))
end 

function add_offset!(self::SWC, offset::Tuple)
    @assert length(offset) == 3
    for in in length( self.points )
        self.points[i].x += offset[1]
        self.points[i].y += offset[2]
        self.points[i].z += offset[3]
    end 
end 

end # module of SWCs

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

function Base.String(self::PointObj)
    "$(self.point_type) $(self.x) $(self.y) $(self.z) $(self.radius) $(self.parent)"
end 

typealias SWC Vector{PointObj}

function SWC(nodeNet::NodeNet)
    edges = NodeNets.get_edges(nodeNet)
    swc = SWC()
    sizehint!(swc, NodeNets.get_node_num(nodeNet))

    for node in NodeNets.get_nodes(nodeNet)
        point = PointObj(0, node[1], node[2], node[3], node[4], -1)
        push!(swc, point)
    end
    # assign parents according to edge 
    for e in edges 
        swc[e[2]].parent = e[1]
    end  
    swc
end 

################## IO ###############################
function save(self::SWC, file_name::AbstractString)
    f = open(file_name, "w")
    for i in 1:length(self)
        write(f, "$i $(String(swc[i])) \n")
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
note that only stretch the coordinates here, not including the radius
since radius was already adjusted in the neighborhood weights 
"""
function stretch_coordinates!(self::SWC, expansion::Tuple)
    @assert length(expansion) == 3
    for i in 1:length( self )
        self[i].x       *= expansion[1]
        self[i].y       *= expansion[2]
        self[i].z       *= expansion[3]
        self[i].radius  *= (prod(expansion))^(1/3)
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
    for in in length( self )
        self[i].x += offset[1]
        self[i].y += offset[2]
        self[i].z += offset[3]
    end 
end 

end # module of SWCs

module SWCs
include("PointObjs.jl")
using .PointObjs 
const ONE_UINT32 = UInt32(1)
export SWC

const SWC = Vector{PointObj}


function SWC( swcString::AbstractString )
    swc = SWC()
    for line in split(swcString, "\n")
        if isempty(line)
            continue 
        end 
        try 
            numbers = map(Meta.parse, split(line))
            # construct a point object
            pointObj = PointObj( numbers[2:7]... )
            push!(swc, pointObj)
        catch err 
            if occursin(r"^\s*#", line)
                println("comment in swc file: $line")
            elseif isempty(line) || line == " " || line == "\n"
                continue 
            else
                @warn("invalid line: $line")
            end 
        end 
    end
    swc
end

################## properties #######################
function get_node_num(self::SWC) length(self) end

function get_edge_num(self::SWC)
    num_edges = 0
    for pointObj in self 
        if PointObjs.get_parent( pointObj ) != -1
            num_edges += 1
        end 
    end 
    num_edges 
end

"""
    get_edges(self::SWC)

get the edges represented as a Vector{NTuple{2,UInt32}}
"""
function get_edges(self::SWC)
    edges = Vector{NTuple{2,UInt32}}()
    for (index, pointObj) in enumerate(self)
        if PointObjs.get_parent( pointObj ) != -1
            push!(edges, ( PointObjs.get_parent( pointObj ), index ))
        end 
    end 
    edges 
end 

"""
    get_total_length( self::SWC )
accumulate all the euclidean distance of edges 
"""
function get_total_length( self::SWC )
    total_length = Float64(0)
    edges = get_edges( self )
    for e in edges 
        total_length += PointObjs.euclidean_distance( self[e[1]], self[e[2]] )
    end 
    total_length 
end 

################## IO ###############################

"""
get binary buffer formatted as neuroglancer nodeNet.

# Binary format
    UInt32: number of vertex
    UInt32: number of edges
    Array{Float32,2}: Nx3 array, xyz coordinates of vertex
    Array{UInt32,2}: Mx2 arrray, node index pair of edges
reference: 
https://github.com/seung-lab/neuroglancer/wiki/Skeletons
"""
function get_neuroglancer_precomputed(self::SWC)
    @show get_node_num(self)
    @show get_edge_num(self)
    # total number of bytes
    num_bytes = 4 + 4 + 4*3*get_node_num(self) + 4*2*get_edge_num(self)
    buffer = IOBuffer( read=false, write=true, maxsize=num_bytes )
    # write the number of vertex, and edges
    write(buffer, UInt32(get_node_num(self)))
    write(buffer, UInt32(get_edge_num(self)))
    # write the node coordinates
    for pointObj in self 
        write(buffer, pointObj.x)
        write(buffer, pointObj.y)
        write(buffer, pointObj.z)
    end
    # write the edges
    for edge in get_edges( self )
        # neuroglancer index is 0-based
        write(buffer, UInt32( edge[1]-ONE_UINT32 ))
        write(buffer, UInt32( edge[2]-ONE_UINT32 ))
    end
    bin = Vector{UInt8}(take!(buffer))
    close(buffer)
    return bin 
end 

function Base.String(self::SWC)
    io = IOBuffer()
    for (index, pointObj) in enumerate(self)
        write(io, "$index $(String(pointObj)) \n")
    end 
    str = String(take!(io))
    close(io)
    str 
end 

function deserialize(data::Vector{UInt8})
    # a pointObj is 21 byte
    @assert mod(length(data), 21) == 0 "the binary file do not match the byte layout of pointObj."
    nodeNum = div(length(data), 21) 
    swc = Vector{PointObj}()
    sizehint!(swc, nodeNum)
    for i in 1:nodeNum 
        pointObj = PointObj( data[(i-1)*21+1:i*21] )
        push!(swc, pointObj)
    end 
    swc
end 

"""
    load_swc_bin( fileName::AbstractString )
load binary swc file 
"""
function load_swc_bin( fileName::AbstractString )
    read( fileName ) |> deserialize
end 

function load_swc(fileName::AbstractString)
    swcString = read( fileName , String)
    SWC( swcString )    
end

function load(fileName::AbstractString)
    if endswith(fileName, ".swc")
        return load_swc(fileName)
    elseif endswith(fileName, ".swc.bin")
        return load_swc_bin(fileName)
    else 
        error("only support .swc and .swc.bin, but this file is: ", fileName)
    end 
end 

function serialize(self::SWC)
    io = IOBuffer( read=false, write=true, maxsize=length(self)*21 )
    for pointObj in self
        write(io, PointObjs.serialize(pointObj))
    end
    data = take!(io)   
end 

"""
    save_swc_bin( self::SWC, fileName::AbstractString )
represent swc file as binary file. the data structure is the same with swc.
"""
function save_swc_bin( self::SWC, fileName::AbstractString )
    data = serialize(self)
    write(fileName, data)
end 

function save_swc(self::SWC, file_name::AbstractString)
    f = open(file_name, "w")
    for (index, pointObj) in enumerate(self)
        write(f, "$index $(String(pointObj)) \n")
    end
    close(f)
end 

function save(self::SWC, fileName::AbstractString)
    if endswith(fileName, ".swc")
        save_swc(self, fileName)
    elseif endswith(fileName, ".swc.bin")
        save_swc_bin(self, fileName)
    else 
        error("only support .swc and .swc.bin, but this file is: ", fileName)
    end 
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
    stretch_coordinates!(self, (2^(mip-1), 2^(mip-1), 1))
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

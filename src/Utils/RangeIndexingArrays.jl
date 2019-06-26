module RangeIndexingArrays

using Statistics
using DataFrames 
using CSV

export RangeIndexingArray

"""
RangeIndexingArray 
"""
struct RangeIndexingArray{T, N} 
    # the septa values defining the indexing ranges 
    septaListTuple::NTuple{N}
    # the table to be indexed 
    table::AbstractArray{T, N}
end 


"""
    range_string_list2septa_list(rangeStringList::Vector) 
transform list of string to a list of septas 

Example:
    input:  3-element Array{Union{Missing, String},1}:
            [ "(0,0.75]",  "(0.75,1.5]",  "(1.5,2]" ]
    
    output: 4-element Array{Float32,1}:
            [  0.0 ,    0.75,    1.5,    2.0 ]
"""
function range_string_list2septa_list(rangeStringList::Vector)
    septaList = Vector{Float32}(undef, length(rangeStringList)+1)
    @inbounds for (i, distRangeString) in rangeStringList |> enumerate
        distRangeString = replace(distRangeString, "("=>"")
        distRangeString = replace(distRangeString, "]"=>"")
        start, stop = map(x->Meta.parse(x)|>Float32, split(distRangeString, ","))
        septaList[i] = start
        if i == length(rangeStringList)
            septaList[end] = stop
        end
    end
    # make sure that the septa were sorted assendingly 
    # so we can use binary search algorithm directly 
    @assert all(diff(septaList) .> zero(Float32))
    septaList
end 

@inline function range_string_list2septa_list(rangeStringList::CSV.Column{String, String})
    range_string_list2septa_list(map(string, rangeStringList))
end 


"""
    RangeIndexingArray{T}(df::DataFrame) 

The dataframe is a matrix, so the result is 2D 
"""
function RangeIndexingArray{T}(df::DataFrame) where {T}
    table = convert( Matrix{T}, df[1:end, 2:end] )
    
    septaList1 = range_string_list2septa_list(df[:Column1] |> Vector{String})
    rangeStringList = map(string, names(df)[2:end] )
    septaList2 = range_string_list2septa_list( rangeStringList )
    
    RangeIndexingArray{T, 2}(tuple(septaList1, septaList2), table)
end 

"""
    RangeIndexingArray(; fileName::AbstractString=joinpath(@__DIR__, "../asset/zfish_score_table.csv"))

read default zfish score table
"""
function RangeIndexingArray{T}(; fileName::AbstractString=joinpath(@__DIR__, "../../asset/zfish_score_table.csv")) where {T}
    df = CSV.read(fileName)
    RangeIndexingArray{T}(df)
end 

"""
    DataFrame(ria::RangeIndexingArray{T,N})
transform a RangeIndexingArray to DataFrame 
"""
function DataFrame(ria::RangeIndexingArray{T,N}) where {T,N}
    distRangeList = Vector{Symbol}()
    @inbounds for i in 1:length(ria.septaListTuple[1])-1
        start = ria.septaListTuple[1][i]
        stop = ria.septaListTuple[1][i+1]
        range = Symbol("($(start),$(stop)]")
        push!(distRangeList, range)
    end

    df = DataFrame(Column1=distRangeList)

    @inbounds for i in 1:length(ria.septaListTuple[2])-1
        start = ria.septaListTuple[2][i]
        stop = ria.septaListTuple[2][i+1]
        range = Symbol("($(start),$(stop)]")
        df[range] = ria.table[:,i]
    end

    df
end 

function Base.eltype(self::RangeIndexingArray{T,N}) where {T,N}
    return T
end

@inline function get_septa_list_tuple(self::RangeIndexingArray)
    self.septaListTuple
end 

@inline function get_table_array(self::RangeIndexingArray)
    self.table 
end 

"""
    float_index2table_index(floatIndex::T, septaList::Vector{T})
transform index with type of floating point to the integer index in the table 
To-Do: use the [Binary Search Algorithm](https://en.wikipedia.org/wiki/Binary_search_algorithm)
"""
function _float_index2table_index(floatIndex::T, septaList::Vector{T}) where T
    if floatIndex<septaList[1]
        return 1
    end 
    @inbounds for i in 1:length(septaList)-1
        if floatIndex <= septaList[i+1] 
            return i 
        end 
    end 
    # this is even larger than the last septa value 
    return length(septaList)-1
end 

@inline function Base.fill!(self::RangeIndexingArray{T,N}, x::T) where {T,N}
    fill!(self.table, x)
end 

function Base.getindex(self::RangeIndexingArray{T,N}, idxes::Integer...) where {T,N}
    table = get_table_array(self)
    return table[idxes...]
end 

"""
    Base.getindex(self::RangeIndexingArray{T,N}, idx::T...)

given a index which could be float type, return a value in the data array. 
"""
function Base.getindex(self::RangeIndexingArray{T,N}, floatIndexList::AbstractFloat...) where {T,N}
    tableIndexList = map((x,y)->_float_index2table_index(x,y), 
                         floatIndexList, get_septa_list_tuple(self))
    table = get_table_array(self)

#    for (i,index) in tableIndexList |> enumerate 
#        if index<1
#            # if they are pretty similar, return score of 1
#            return one(T)
#        end  
#    end 
#    for (i,index) in tableIndexList |> enumerate 
#        if index>size(table, i)
#            # if they are very different, return score of 0
#            return zero(T)
#        end 
#    end 

    return get_table_array(self)[tableIndexList...,]
end 

function Base.setindex!(self::RangeIndexingArray{T,N}, value::T, 
                       inds::AbstractFloat...) where {T,N}
    tableIndexList = map((x,y)->_float_index2table_index(x,y), inds, get_septa_list_tuple(self))
    table = get_table_array(self)
    table[tableIndexList...,] = value 
end

"""
    to_position_only(self::RangeIndexingArray{T,2})

convert to position only range indexing array. This is used for dendrite query. 
The dendrite seems only respect approximity rather than direction.
"""
function to_position_only(self::RangeIndexingArray{T,2}) where T
    ret = deepcopy(self)
    v = mean(ret.table; dims=2) |> vec
    for j in size(ret.table, 2)
        ret.table[:,j] = v
    end
    @assert size(ret.table, 2) == 10
    ret
end

end # end of module 

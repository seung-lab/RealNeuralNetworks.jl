#!/usr/bin/env julia

module DBFs

const DBF = Vector{Float32}

export DBF 

using Base.Cartesian

"""
use segmentation to get binary image to save memory usage
"""
function compute_DBF(seg::Array{T,3}, obj_id::T) where T
    error("unimplemented")
end 

"""

    compute_DBF( pointCloud )

  Returns an array of DBF values for the point cloud. Currently creates
  a binary image, and runs bwd2 on it, though ideally we'd get rid of the
  need for an explicit bin_im
"""
function compute_DBF( pointCloud::Array{T, 2} ) where T
    bin_im = create_binary_image( pointCloud );
    compute_DBF(pointCloud, bin_im)
end

"""
    compute_DBF( bin_im )
"""
function compute_DBF( pointCloud::Array{T,2}, bin_im::Union{BitArray, Array{Bool, 3}} ) where T
    dbf_im = distance_transform( bin_im );
    return extract_dbf_values( dbf_im, pointCloud );
end 

"""
compute Distance from Boundary Field (DBF) based on point cloud and the boundary points

WARN: this function do not work correctly!
"""
function compute_DBF( points::Array{T,2}, boundary_point_indexes::Vector ) where T
    error("this function do not work correctly, have a bug!")
    num = size(points, 1)
    dbf = DBF(num)
    fill!(dbf, Inf32)
    for i in 1:num
        point = points[i,:]
        for bpi in boundary_point_indexes
            boundary = points[bpi, :]
            # filter out some far away boundary points 
            if  abs(point[1]-boundary[1]) < MAX_BOUNDARY_DISTANCE && 
                abs(point[2]-boundary[2]) < MAX_BOUNDARY_DISTANCE && 
                abs(point[3]-boundary[3]) < MAX_BOUNDARY_DISTANCE 
                # compute euclidean distance
                ed = norm(point .- boundary)
                if ed < dbf[i]
                    dbf[i] = ed 
                end 
            end 
        end 
    end
    return dbf
end 


#NOTE voxelSize not currently functional
"""

    distance_transform( d::AbstractArray{T,N}, voxelSize::Vector{Float32}=ones(Float32, N) )

  Returns a euclidean distance transformation of the mask provided by d. The return
  value will be a volume of the same size as d where the value at each index corresponds
  to the distance between that location and the nearest location for which d > 0.
"""
@generated function distance_transform( d::AbstractArray{T,N}, voxelSize::Vector{Float32}=ones(Float32, N) ) where {T,N}
  quote

  @assert length(voxelSize) == $N;

  res = zeros(Float32,size(d));

  fill_f0!(res, d);

  dims = 1:$N;

  for d in dims
    vol_voronoi_edt!( res, voxelSize[d] );
    res = permutedims( res, circshift(dims,1) );
  end

  sqrt!(res)
  res

  end#quote
end


"""
Fills an n-dimensional volume with initial states for edt transformation,
inf for non-feature voxels, and 0 for feature voxels
"""
@generated function fill_f0!( arr::Array{Float32,N}, fv::AbstractArray{T,N} ) where {T,N}
  quote

  #apparently this generates lots of allocations,
  # I'm still not entirely sure why... but cool!
  # (learned the answer - choice between 0 and Inf is type unstable)
  @nloops $N i arr begin
    (@nref $N arr i) = (@nref $N fv i) > 0 ? 0. : Inf;
  end

  #arr[fv] = 0;
  #arr[!fv] = Inf;

  end#quote
end


"""
Performs the edt transformation along the first dimension of the N-dimensional
volume
"""
@generated function vol_voronoi_edt!( arr::Array{Float32,N}, dim::Float32 ) where N
  quote

  s1 = size(arr,1);
  g = zeros(Float32,(s1,));
  h = zeros(Int,    (s1,));
  @nloops $N i j->(j==1 ? 0 : 1:size(arr,j)) begin

      fill!(h, zero(eltype(h))); 
      fill!(g, zero(eltype(g)));
    later_indices = (@ntuple $N j->i_j);
    row_voronoi_edt!( arr, later_indices[2:end], g,h, dim ) #F_{d-1}

  end

  end#quote
end


"""
Performs the edt over a specific row in the volume, following the first dimension
"""
@generated function row_voronoi_edt!( F::Array{Float32,N}, indices::Tuple,
  g::Vector{Float32}, h::Vector{Int}, dim::Float32 ) where N
  quote

  #count of potential feature vectors
  numPotentialFeatureVectors::Int = 0;

  #selecting out the value in the row
  @inbounds f = @nref $N F j->(j==1 ? (:) : indices[j-1]);

  #construct set of feature voxels whose voronoi
  # cells intersect the row
  for i in eachindex( f )

    #scanning for possible feature vector locations
    if !isinf( f[i] )
      if numPotentialFeatureVectors < 2
        #l += 1; g[l] = f[i]; h[l] = i;
        numPotentialFeatureVectors += 1; 
        g[numPotentialFeatureVectors] = f[i]; 
        h[numPotentialFeatureVectors] = i;
      else
        while (numPotentialFeatureVectors >= 2 && 
               remove_euclidean_distance_transform(g[numPotentialFeatureVectors-1],
                                                   g[numPotentialFeatureVectors], f[i],
                                                   h[numPotentialFeatureVectors-1],
                                                   h[numPotentialFeatureVectors], i))
          numPotentialFeatureVectors -= 1;
        end
        numPotentialFeatureVectors += 1; 
        g[numPotentialFeatureVectors] = f[i]; 
        h[numPotentialFeatureVectors] = i;
      end #if numPotentialFeatureVectors
    end #if !isinf

  end #for i

  # if no possible feature voxels, stop now
  if numPotentialFeatureVectors == 0; return nothing end

  #assign new closest feature vectors
  # and update new distances
  num_fvs = numPotentialFeatureVectors; numPotentialFeatureVectors = 1;
  for i in eachindex( f )

    #if we haven't reached the end, and the next feature vector is closer
    # to this location than the current one
    while (numPotentialFeatureVectors < num_fvs && 
           (g[numPotentialFeatureVectors] + (h[numPotentialFeatureVectors] - i)^2 > 
                g[numPotentialFeatureVectors+1] + (h[numPotentialFeatureVectors+1] - i)^2))
      numPotentialFeatureVectors += 1;
    end

    (@nref $N F j->(j==1 ? i : indices[j-1])) =  g[numPotentialFeatureVectors] + 
                                            (h[numPotentialFeatureVectors] - i)^2;
  end #for i

  end #quote
end

"""
Getting too tired to document these next few, but will be worth it if it works
"""
function fv_isfurther(g1::Float32 ,h1::Int, g2::Float32, h2::Int, i::Int)
  return g1 + (h1 - i)^2 > g2 + (h2 - i)^2;
end

"""
    remove_euclidean_distance_transform 
"""
function remove_euclidean_distance_transform( g1::Float32, g2::Float32, g3::Float32, h1::Int, h2::Int, h3::Int )
  a = h2 - h1;
  b = h3 - h2;
  c = h3 - h1;

  return (c*g2 - b*g1 - a*g3 - a*b*c > 0)
end


function sqrt!( d::AbstractArray )
  for i in eachindex(d)
    d[i] = sqrt(d[i]);
  end
end


"""

    create_binary_image( pointCloud )

  Creates a boolean volume where the non-segment indices
  map to true, while the segment indices map to false.
"""
function create_binary_image( pointCloud::Array{T,2} ) where T;

  max_dims = maximum( pointCloud, dims=1 );
    bin_im = falses(max_dims...)

  for p in 1:size( pointCloud, 1 )
    bin_im[ pointCloud[p,:]... ] = false;
  end

  bin_im;
end

"""
    create_binary_image( seg, obj_id )

Creates a boolean volume where the non-segment indices
map to true, while the segment indices map to false 
"""
function create_binary_image( seg::Array{T,3}; obj_id::T = one(T) ) where T
    bin_im = trues(size(seg))
    for i in eachindex(seg)
        if seg[i] == obj_id 
            bin_im[i] = false 
        end 
    end 
    bin_im
end 


"""

    extract_dbf_values( dbf_image, pointCloud )

  Takes an array where rows indicate subscripts, and extracts the values
  within a volume at those subscripts (in row order)
"""
@generated function extract_dbf_values( dbf_image::Array{Float32,N}, pointCloud ) where N
  quote

  num_points = size( pointCloud, 1 );
  dbf_values = zeros(Float32, num_points);

  for p in 1:num_points 
    dbf_values[p] = (@nref $N dbf_image i->pointCloud[p,i]);
  end

  dbf_values
  end#quote
end


end#module

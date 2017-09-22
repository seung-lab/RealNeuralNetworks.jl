#!/usr/bin/env julia

__precompile__()

module DBFs

using Base.Cartesian

typealias DBF Vector{Float32}

"""
use segmentation to get binary image to save memory usage
"""
function compute_DBF{T}(point_cloud::Array{UInt32,2}, seg::Array{UInt32,3}, obj_id::T)
    bin_im = create_boundary_image( seg, obj_id ) 
    compute_DBF( point_cloud, bin_im )
end 

"""

    compute_DBF( point_cloud )

  Returns an array of DBF values for the point cloud. Currently creates
  a binary image, and runs bwd2 on it, though ideally we'd get rid of the
  need for an explicit bin_im
"""
function compute_DBF{T}( point_cloud::Array{T, 2} )
    bin_im = create_boundary_image( point_cloud );
    compute_DBF(point_cloud, bin_im)
end

"""
    compute_DBF( bin_im )
"""
function compute_DBF{T}( point_cloud::Array{T,2}, bin_im::Array{Bool, 3} )
    dbf_im = bwd2( bin_im );
    dbf_vec = extract_dbf_values( dbf_im, point_cloud );
    dbf_vec;
end 

"""
compute Distance from Boundary Field (DBF) based on point cloud and the boundary points

WARN: this function do not work correctly!
"""
function compute_DBF{T}( points::Array{T,2}, boundary_point_indexes::Vector )
    error("this function do not work correctly, have a bug!")
    num = size(points, 1)
    dbf = Vector{Float32}(num)
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


#NOTE voxel_dims not currently functional
"""

    bwd2( d::AbstractArray{T,N}, voxel_dims::Vector{Float64}=ones(N) )

  Returns a euclidean distance transformation of the mask provided by d. The return
  value will be a volume of the same size as d where the value at each index corresponds
  to the distance between that location and the nearest location for which d > 0.
"""
@generated function bwd2{T,N}( d::AbstractArray{T,N},
  voxel_dims::Vector{Float64}=ones(N) )
  quote

  @assert length(voxel_dims) == $N;

  res = zeros(Float64,size(d));

  fill_f0!(res, d);

  dims = 1:$N;

  for d in dims
    vol_voronoi_edt!( res, voxel_dims[d] );
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
@generated function fill_f0!{T,N}( arr::Array{Float64,N}, fv::AbstractArray{T,N} )
  quote

  #apparently this generates lots of allocations,
  # I'm still not entirely sure why... but cool!
  # (learned the answer - choice between 0 and Inf is type unstable)
  @nloops $N i arr begin
    (@nref $N arr i) = (@nref $N fv i) > 0? 0. : Inf;
  end

  #arr[fv] = 0;
  #arr[!fv] = Inf;

  end#quote
end


"""
Performs the edt transformation along the first dimension of the N-dimensional
volume
"""
@generated function vol_voronoi_edt!{N}( arr::Array{Float64,N}, dim::Float64 )
  quote

  s1 = size(arr,1);
  g = zeros(Float64,(s1,));
  h = zeros(Int,    (s1,));
  @nloops $N i j->(j==1?0:1:size(arr,j)) begin

    fill!(h,0); fill!(g,0);
    later_indices = (@ntuple $N j->i_j);
    row_voronoi_edt!( arr, later_indices[2:end], g,h, dim ) #F_{d-1}

  end

  end#quote
end


"""
Performs the edt over a specific row in the volume, following the first dimension
"""
@generated function row_voronoi_edt!{N}( F::Array{Float64,N}, indices::Tuple,
  g::Vector{Float64}, h::Vector{Int}, dim::Float64 )
  quote

  #count of potential feature vectors
  l::Int = 0;

  #selecting out the value in the row
  @inbounds f = @nref $N F j->(j==1?(:):indices[j-1]);

  #construct set of feature voxels whose voronoi
  # cells intersect the row
  for i in eachindex( f )

    #scanning for possible feature vector locations
    if !isinf( f[i] )
      if l < 2
        #l += 1; g[l] = f[i]; h[l] = i;
        l += 1; g[l] = f[i]; h[l] = i;
      else
        while (l >= 2 && removeEDT(g[l-1],g[l],f[i],h[l-1],h[l],i))
          l -= 1;
        end
        l += 1; g[l] = f[i]; h[l] = i;
      end #if l
    end #if !isinf

  end #for i

  # if no possible feature voxels, stop now
  if l == 0; return nothing end

  #assign new closest feature vectors
  # and update new distances
  num_fvs = l; l = 1;
  for i in eachindex( f )

    #if we haven't reached the end, and the next feature vector is closer
    # to this location than the current one
    while (l < num_fvs && (g[l] + (h[l] - i)^2 > g[l+1] + (h[l+1] - i)^2))
      l += 1;
    end

    (@nref $N F j->(j==1?i:indices[j-1])) = g[l] + (h[l] - i)^2;
  end #for i

  end #quote
end

"""
Getting too tired to document these next few, but will be worth it if it works
"""
function fv_isfurther(g1::Float64 ,h1::Int, g2::Float64, h2::Int, i::Int)
  return g1 + (h1 - i)^2 > g2 + (h2 - i)^2;
end

"""
    removeEDT
"""
function removeEDT( g1::Float64, g2::Float64, g3::Float64, h1::Int, h2::Int, h3::Int )
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

    create_boundary_image( point_cloud )

  Creates a boolean volume where the non-segment indices
  map to true, while the segment indices map to false.
"""
function create_boundary_image{T}( point_cloud::Array{T,2} );

  max_dims = maximum( point_cloud, 1 );

  bin_im = ones(Bool, max_dims...);

  for p in 1:size( point_cloud, 1 )
    bin_im[ point_cloud[p,:]... ] = false;
  end

  bin_im;
end

"""
    create_boundary_image( seg, obj_id )

Creates a boolean volume where the non-segment indices
map to true, while the segment indices map to false 
"""
function create_boundary_image{T}( seg::Array{T,3}, obj_id::T )
    bin_im = ones(Bool, size(seg))
    for i in eachindex(seg)
        if seg[i] == obj_id 
            bin_im[i] = false 
        end 
    end 
    bin_im
end 


"""

    extract_dbf_values( dbf_image, point_cloud )

  Takes an array where rows indicate subscripts, and extracts the values
  within a volume at those subscripts (in row order)
"""
@generated function extract_dbf_values{N}( dbf_image::Array{Float64,N}, point_cloud )
  quote

  num_points = size( point_cloud, 1 );
  dbf_values = zeros(num_points);

  for p in 1:size(point_cloud, 1)
    dbf_values[p] = (@nref $N dbf_image i->point_cloud[p,i]);
  end

  dbf_values
  end#quote
end


end#module

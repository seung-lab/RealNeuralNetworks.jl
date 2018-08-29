module FakeSegmentations
using LinearAlgebra

const ZERO_UINT32 = UInt32(0)
const ONE_UINT32 = UInt32(1)
const RADIUS = Float32(10)
const HEIGHT = 100
const THICKNESS = 10

function cylinder(;radius::Float32=RADIUS, height::Integer=HEIGHT)
    width = Int( radius*2 + 1 )
    center = Int(radius + 1)
    seg = zeros(UInt32, (width,width,height))
    for x in Int(center-radius):Int(center+radius)
        for y in Int(center-radius):Int(center+radius)
            if norm([x,y].-[center, center]) < radius 
                seg[x,y,:] .= ONE_UINT32 
            end 
        end 
    end 
    seg
end

function broken_cylinder(;radius::Float32=RADIUS, height::Integer=HEIGHT, 
                         breakStartZ::Integer=div(height,2), 
                         thickness::Integer = THICKNESS )
    seg = cylinder(radius=radius, height=height)
    # create a breaking part
    seg[ :,:, breakStartZ : breakStartZ+10 ] .= ZERO_UINT32
    seg
end 

function ring(; centerLineRadius::Integer = 50, ringRadius::Integer=5)
    seg = zeros(UInt32, (2*centerLineRadius+2*ringRadius + 5, 
                         2*centerLineRadius+2*ringRadius + 5, 2*ringRadius + 5) )
    center = [map(x->div(x,2), size(seg))...]
    for y in 1:size(seg, 2)
        for x in 1:size(seg, 1)
            d = norm([x,y] .- center[1:2])
            if abs(d-centerLineRadius) < 2 
                # this point is on the center line of the ring
                for rz in center[3]-ringRadius : center[3]+ringRadius 
                    for ry in y-ringRadius : y+ringRadius 
                        for rx in x-ringRadius : x+ringRadius
                            d2 = norm([x,y,center[3]] .- [rx,ry,rz])
                            if d2 <= ringRadius
                                seg[rx,ry,rz] = ONE_UINT32 
                            end 
                        end 
                    end 
                end 
            end 
        end 
    end 
    @assert ndims(seg) == 3
    seg
end 

function broken_ring(; centerLineRadius::Integer = 50, ringRadius::Integer=5)
    seg = ring(; centerLineRadius = centerLineRadius, ringRadius = ringRadius)
    center = map(x->div(x,2), size(seg))
    d = div(centerLineRadius, 10)
    seg[center[1]-d:center[1]+d, :, :] .= ZERO_UINT32 
    seg[:, center[2]-d:center[2]+d, :] .= ZERO_UINT32 
    seg
end 

function Y_shape()
    error("unimplemented")
end 

end # module 

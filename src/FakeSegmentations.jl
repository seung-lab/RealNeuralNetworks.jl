module FakeSegmentations

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
                seg[x,y,:] = ONE_UINT32 
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
    seg[ :,:, breakStartZ : breakStartZ+10 ] = ZERO_UINT32
    seg
end 

function ring(; centerLineRadius::Float32 = Float32(50), ringRadius::Float32=Float32(5))
    error("unimplemented")
end 

function Y_shape()
    error("unimplemented")
end 

end # module 

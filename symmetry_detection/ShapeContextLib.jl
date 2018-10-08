module ShapeContextLib

using LinearAlgebra
using FFTW

export ShapeContext

function XYZToZRT(point, ex, ey, ez)
    x = point[1]*ex[1] + point[2]*ex[2] + point[3]*ex[3]
    y = point[1]*ey[1] + point[2]*ey[2] + point[3]*ey[3]
    z = point[1]*ez[1] + point[2]*ez[2] + point[3]*ez[3]
    z = abs(z)
    theta = atan(y, x)
    if theta < 0.0 theta += 2*pi end
    radius = sqrt(x*x + y*y)
    return z, radius, theta
end

function discretize!(dest, x, bin_size)
    # Discretize a 1d array elementwise by floor(x, bin_size) + 1.
    for i = 1:size(x, 1)
        dest[i] = floor(Int32, x[i] / bin_size) + 1
    end
    return nothing
end

#shape context descriptor
function ShapeContext(points, radius=0.6, rbin=5, hbin=5, abin=128)
    Descriptor = zeros(size(points, 1), rbin*hbin*abin)
    points_temp = zeros(size(points, 1),3)
    r = zeros(size(points, 1),1)
    height = zeros(size(points, 1),1)
    radi = zeros(size(points, 1),1)
    angle = zeros(size(points, 1),1)
    height_q = zeros(Int32, size(points, 1),1)
    radi_q = zeros(Int32, size(points, 1),1)
    angle_q = zeros(Int32, size(points, 1),1)
    hist = zeros(abin, hbin, rbin)
    Fhist = zeros(ComplexF32, abin, hbin, rbin)
    for p = 1:size(points, 1)
        broadcast!(-, points_temp, points, view(points, p:p, :))

        ez = points[p,:] / norm(points[p,:])
        ex = cross(Float64[0,0,1], vec(ez))
        if norm(ex)<0.001
            ex = Float64[1,0,0]
        else
            ex = ex / norm(ex)
        end
        ey = cross(vec(ez), ex)
        
        for i = 1:size(points, 1)
            r[i] = sqrt(points_temp[i,1]^2+points_temp[i,2]^2+points_temp[i,3]^2)
            height[i], radi[i], angle[i] = XYZToZRT(view(points_temp, i, :), ex, ey, ez)
        end

        discretize!(height_q, height,  radius/hbin)
        discretize!(radi_q, radi, radius/rbin)
        discretize!(angle_q, angle,  2*pi/abin)

        for i in eachindex(hist) hist[i]=0 end
        count = 0
        for i = 1:size(points_temp, 1)
            if r[i] > 0.02 && r[i] < radius && height_q[i] <= hbin && angle_q[i] <= abin && radi_q[i] <= rbin
                hist[angle_q[i], height_q[i], radi_q[i]] += 1
                count += 1
            end
        end
        for i in eachindex(hist)
            Fhist[i] = complex(hist[i] / count)
        end
        fft!(Fhist, 1)
        for i in eachindex(Fhist)
            hist[i] = abs(Fhist[i])
        end
        setindex!(Descriptor, hist, p, :)  # fast way of doing Descriptor[p,:]=hist[:]
    end
    Descriptor
end


end # module end


using MeshIO
using FileIO
using Distances
using NearestNeighbors
using ICPUtil
@everywhere include("utils.jl")

function XYZToZRT(point, ex, ey, ez)
    x = point[1]*ex[1] + point[2]*ex[2] + point[3]*ex[3]
    y = point[1]*ey[1] + point[2]*ey[2] + point[3]*ey[3]
    z = point[1]*ez[1] + point[2]*ez[2] + point[3]*ez[3]
    z = abs(z)
    theta = atan2(x, y)
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
@everywhere function ShapeContext(points, radius=0.6, rbin=5, hbin=5, abin=128)
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
    Fhist = zeros(Complex64, abin, hbin, rbin)
    for p = 1:size(points, 1)
        broadcast!(-, points_temp, points, sub(points, p, :))

        ez = points[p,:] ./ norm(points[p,:])
        ex = cross(Float64[0,0,1], vec(ez))
        if norm(ex)<0.001
            ex = Float64[1,0,0]
        else
            ex = ex ./ norm(ex)
        end
        ey = cross(vec(ez), ex)
        
        for i = 1:size(points, 1)
            r[i] = sqrt(points_temp[i,1]^2+points_temp[i,2]^2+points_temp[i,3]^2)
            height[i], radi[i], angle[i] = XYZToZRT(sub(points_temp, i, :), ex, ey, ez)
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

@everywhere function calculateTransformation(reflection, point_x1, point_x2, point_y1, point_y2)
    point_x = point_x1 ./ norm(point_x1)
    point_y = point_y1 ./ norm(point_y1)
    diff_x = point_x2 - point_x .* (point_x2 * point_x')
    diff_y = point_y2 - point_y .* (point_y2 * point_y')
    diff_x = diff_x ./ norm(diff_x)
    diff_y = diff_y ./ norm(diff_y)
    cross_x = cross(vec(point_x), vec(diff_x))
    cross_y = cross(vec(point_y), vec(diff_y))
    if reflection==1
        cross_y = -cross_y
    end
    transy = vcat(point_y, diff_y, cross_y')
    transx =  vcat(point_x, diff_x, cross_x')
    transy' * transx
end 

@everywhere function processSymmetry(symtrans, symscores)
    axises = []
    angles = []
    reflections = []
    for i = 1:size(symtrans, 1)
        axis, ang, reflection = matrix2axis(symtrans[i])
        push!(axises, axis)
        push!(angles, ang)
        push!(reflections, reflection)
    end
    idx = sortperm(symscores)
    axises = axises[idx]
    angles = angles[idx]
    reflections = reflections[idx]
    axises_f = []
    angles_f = Array(Any,0)
    reflection_f = []
    inversion_f = 1
    axisnorm = [norm(axises[i]) for i = 1:size(axises,1)]
    idx = find(axisnorm .< 0.1)
    inversionidx = find(reflections[idx]==-1)
    if size(inversionidx,1) > 0 inversion_f = -1 end
    idx = axisnorm .> 0.1
    axises = axises[idx]
    angles = angles[idx]
    reflections = reflections[idx]
    for i = 1:size(axises, 1)
        if reflections[i] == 1 continue end
        for j = i+1:size(axises, 1)
            if reflections[j] == 1 continue end
            axis = cross(vec(axises[i]), vec(axises[j]))
            push!(axises, axis ./ norm(axis))
            push!(angles, acos(dot(vec(axises[i]), vec(axises[j])))*2)
            push!(reflections, 1)
            #println(size(axises, 1))
            #println(size(angles, 1))
        end
    end
    if size(axises,1)>0
        axises_f = axises[1]
        push!(angles_f,Float64[angles[1]])
        reflection_f = Int64[reflections[1]]
        for i = 2:size(axises,1)
            dist = abs(axises[i]'*axises_f)
            maxval, maxidx = findmax(dist)
            if maxval > 0.97
                if reflections[i]==1
                    push!(angles_f[maxidx],angles[i])
                    elseif angles[i] < 0.01
                    reflection_f[maxidx]= -1
                end
            else
                axises_f = hcat(axises_f, axises[i])
                if reflections[i]==1
                    push!(angles_f,Float64[angles[i]])
                    reflection_f = hcat(reflection_f, 1)
                elseif angles[i] < 0.01
                    push!(angles_f,Float64[0])
                    reflection_f = hcat(reflection_f, -1)
                else
                    push!(angles_f,Float64[0])
                    reflection_f = hcat(reflection_f, 1)
                end
            end
        end
    end

    #print(axises_f)
    #print(angles_f)
    #print(reflection_f)

    angles_gt = Array(Any,0)
    push!(angles_gt, [0])
    push!(angles_gt, [0, pi])
    push!(angles_gt, [0, 2/3*pi])
    push!(angles_gt, [0, pi/2, pi])
    push!(angles_gt, [0, 2/5*pi, 4/5*pi])
    push!(angles_gt, [0, pi/3, 2*pi/3, pi])
    push!(angles_gt, [0, 2/7*pi, 4/7*pi, 6/7*pi])
    push!(angles_gt, [0, pi/4, pi/2, pi*3/4, pi])

    if isempty(axises_f)
        return ([],[],[])
    else
        degree_f = zeros(Int, 1, size(axises_f,2))
        for i = 1:size(axises_f,2)
            for d = 1:8
                diff = pairwise(Cityblock(), angles_f[i]', angles_gt[d]')
                dist = maximum(minimum(diff,2),1)
                if dist[1] < 0.05
                    degree_f[i]=d
                    break
                end
            end
        end
    end

    axises_f, degree_f, reflection_f
end

@everywhere function saveSymmetryAxis(outName::AbstractString, axises_f, degree_f, reflection_f)
    fout = open(outName, "w")
    for i = 1:3
        for j = 1:size(axises_f,2)
            @printf(fout, "%7.3f ", axises_f[i,j])
        end
        @printf(fout, "\n")
    end
    for j = 1:size(degree_f,2)
        @printf(fout, "%3d ", degree_f[j])
    end
    @printf(fout, "\n")
    for j = 1:size(reflection_f,2)
        @printf(fout, "%3d ", reflection_f[j])
    end
    @printf(fout, "\n")
    close(fout)
    return 0
end

@everywhere function refineTransform(candtrans, points)
    symtrans = []
    symscores = []
    for transform in candtrans
        points_trans = transform * points'
    #    TR,TT,ER,t = icp(points', points_trans, 5)
    #    transform = TR * transform
    #    if ER[end] < 0.02 
    #        push!(symtrans, transform)
    #        push!(symscores, ER[end])
    #    else 
            TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point")
            transform = TR * transform
            if ER[end] < 0.025
                push!(symtrans, transform)
                push!(symscores, ER[end])
            end
    #    end
    end

    symtrans, symscores
end

@everywhere function symmetrySpectral(points::Array{Float64}, desc::Array{Float64})
    if isempty(desc) || size(desc,1) != size(points,1)
        desc = ShapeContext(points, 1.5, 6,6,128)
    end
    dists = pairwise(Euclidean(), desc')
    S = exp(-dists.^2/0.02)
    println("compute feature done.")
    tic()
    orbitsize = sum(S,2)
    #minsize = minimum(orbitsize)
    #minsize = max(10, minsize)
    numEffectCal = 0
    numMaxCal = 400
    symtrans = []
    symscore = []
    kdtree = KDTree(points')
    pointdists = pairwise(Euclidean(), points')
    validpairs = (dists .< 0.15) & (pointdists .> 0.1)
    istrue(x) = x==true
    pairidx = find(istrue, validpairs)
    if size(pairidx,1) < 200
        println("too little matching pairs, stop");
        return 0;
    end
    candtrans = []
    while (numEffectCal < numMaxCal)
        idx_pair1 = rand(pairidx)
        idx_pair2 = rand(pairidx)
        idx_x1 = floor(Int,(idx_pair1-1)/size(points,1))+1
        idx_y1 = (idx_pair1-1) % size(points,1)+1
        idx_x2 = floor(Int,(idx_pair2-1)/size(points,1))+1
        idx_y2 = (idx_pair2-1) % size(points,1)+1

        if norm(points[idx_x1,:]-points[idx_x2,:]) < 0.05 
            continue; 
        end
        if abs(norm(points[idx_x1,:]-points[idx_x2,:]) - norm(points[idx_y1,:]-points[idx_y2,:])) > 0.1 
            continue; 
        end
        #this two conditions are only for global symmetry
        if abs(norm(points[idx_x1,:])-norm(points[idx_y1,:])) > 0.1 
            continue; 
        end
        if abs(norm(points[idx_x2,:])-norm(points[idx_y2,:])) > 0.1 
            continue; 
        end
        for reflection = 0:1
            trans = calculateTransformation(reflection, points[idx_x1,:], points[idx_x2,:], points[idx_y1,:], points[idx_y2,:])
            points_trans = trans * points'
            match, mindist =knn(kdtree, points_trans, 1)
            if mean(mindist)[1]<0.2
                push!(candtrans, trans)
            end
            numEffectCal +=1
        end
    end
    toc()

    candtrans
end

@everywhere function symmetrySpectral(modelID::AbstractString)
    #generate pointcloud
    pointcloud = load(modelID*".off")
    points = zeros(size(pointcloud.vertices, 1), 3)
    for i in 1:size(pointcloud.vertices,1)
        points[i, :] = [pointcloud.vertices[i][1], pointcloud.vertices[i][2], pointcloud.vertices[i][3]]
    end
    println("read pointcloud")
    @printf("%d points\n", size(points, 1))
    candtrans = symmetrySpectral(points)
    symtrans, symscores = refineTransform(candtrans, points)
    axises_f, degree_f, reflection_f = processSymmetry(symtrans, symscores)

    outName = modelID*".sym"
    saveSymmetryAxis(outName, axises_f, degree_f, reflection_f)
end

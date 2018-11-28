module SymmetrySpectralLib

using Distances
using NearestNeighbors
using ICPUtil
using TransformUtil
using ShapeContextLib

export refineTransform, processSymmetry, symmetrySpectral

function calculateTransformation(reflection, point_x1, point_x2, point_y1, point_y2)
    point_x = point_x1 / norm(point_x1)
    point_y = point_y1 / norm(point_y1)
    diff_x = point_x2 - point_x .* (point_x2 * point_x')
    diff_y = point_y2 - point_y .* (point_y2 * point_y')
    diff_x = diff_x / norm(diff_x)
    diff_y = diff_y / norm(diff_y)
    cross_x = cross(vec(point_x), vec(diff_x))
    cross_y = cross(vec(point_y), vec(diff_y))
    if reflection==1
        cross_y = -cross_y
    end
    transy = vcat(point_y, diff_y, cross_y')
    transx =  vcat(point_x, diff_x, cross_x')
    transy' * transx
end 

function refineTransform(candtrans, points)
    symtrans = []
    symscores = []
    for transform in candtrans
        points_trans = transform * points'
        TR,TT,ER,HD = icp(points', points_trans, 50, Minimize="point")
        transform = TR * transform
        if ER[end] < 0.035
            push!(symtrans, transform)
            push!(symscores, ER[end])
        end
    end

    symtrans, symscores
end

function processSymmetry(symtrans, symscores)
    axises = []
    angles = []
    reflections = []
    scores = []
    for i = 1:size(symtrans, 1)
        axis, ang, reflection = matrix2axis(symtrans[i])
        if ang > 0.15 || reflection == -1
            push!(axises, axis)
            push!(angles, ang)
            push!(reflections, reflection)
            push!(scores, symscores[i])
        end
    end
    idx = sortperm(scores)
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
            push!(axises, axis / norm(axis))
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
            if maxval > 0.90
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

function symmetrySpectral(points::Array{Float64}, desc::Array{Float64})
    if isempty(desc) || size(desc,1) != size(points,1)
        desc = ShapeContext(points, 1.5, 6,6,128)
    end
    dists = pairwise(Euclidean(), desc')
    println("compute feature done.")
    tic()
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


end # module end

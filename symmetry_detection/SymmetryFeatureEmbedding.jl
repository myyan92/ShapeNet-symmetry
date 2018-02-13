numCore = 10
addprocs(numCore - 1)


push!(LOAD_PATH, "../../geometryprocessing/julia")
@everywhere using MeshIO
@everywhere using FileIO
@everywhere using MyGeomUtils
include("./ConnectedComponent.jl")
include("./symmetrySpectral.jl")
include("./icp.jl")

@everywhere function thresh_start()
    return 0.05
end

@everywhere function thresh_end()
    return 0.035
end

#axis is the normal of the reflection plane. 
#return refined axis or none if not a reflection plane.
@everywhere function refineAxis_reflect(points, axis, log)
    @printf(log, "%s\n", "refine_reflect:")
    transform = axis2matrix(axis, -1)
    points_trans = transform * points'
    TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point")
    transform = TR * transform
    @printf(log, "%f %f\n", ER[1], ER[end])
    if ER[1] < thresh_start() && ER[end] < thresh_end()
        axis, angle, reflect = matrix2axis(transform)
        return axis
    end
    return []
end

#axis is the rotation axis
#reflectPose and axis forms one of the reflection plane
@everywhere function refineAxis_C(points, axis, degree, reflectPose, log)
    @printf(log, "%s\n", "refine_C:")
    transform = axis2matrix(axis, degree)
    baserotation = transform
    basereflection = axis2matrix(cross(axis,reflectPose), -1)
    valid_rotate = zeros(degree,1)
    valid_reflect = zeros(degree,1)
    isreflect = true
    for j = 1:degree
        points_trans = baserotation * points'
        TR,TT,ER,t = icp(points', points_trans, 1, Minimize="point")
        @printf(log, "%f ", ER[1])
        if ER[1] < thresh_start()
            TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point")
            @printf(log, "%f", ER[end])
            if ER[end] < thresh_end() valid_rotate[j]=1; end
        end
        @printf(log, "\n")
        points_trans = basereflection * points'
        TR,TT,ER,t = icp(points', points_trans, 1, Minimize="point")
        @printf(log, "%f ", ER[1])
        if ER[1] < thresh_start()
            TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point")
            @printf(log, "%f", ER[end])
            if ER[end] < thresh_end() valid_reflect[j]=1; end  
        end
        @printf(log, "\n")
        baserotation = transform * baserotation
        basereflection = transform * basereflection
    end
    degree_f = 1
    idx = 0
    for j = 1:degree
        if mod(degree,j) != 0 continue; end
        numValid = Int(sum(valid_rotate[j:j:end]))
        if j*numValid == degree 
            degree_f = numValid
            valid_reflect = reshape(valid_reflect, j, numValid)
            numValid2= sum(valid_reflect, 2)
            idx = findfirst(numValid2.==numValid)
            if idx == 0
                isreflect= false
            end    
            break
        end
    end
    axis_f = zeros(3,1)
    if degree_f == 1
        axis_f = axis
    else
        transform = axis2matrix(axis, degree_f)
        baserotation = transform
        for j = 1:degree_f
            points_trans = baserotation * points'
            TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point")
            axis_t, degree_t, reflect_t = matrix2axis(TR*baserotation)
            if dot(vec(axis_f), vec(axis_t)) < 0 axis_t = -axis_t; end
            axis_f = axis_f + axis_t
            baserotation = transform * baserotation
        end
        axis_f = axis_f ./ norm(axis_f)
    end
    if isreflect;
        transform = axis2matrix(axis_f, degree)
        basereflection = axis2matrix(cross(vec(axis_f),vec(reflectPose)), -1)
        for j = 1:idx-1
            basereflection = transform * basereflection
        end
        points_trans = basereflection * points'
        TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point")
        axis_t, degree_t, reflect_t = matrix2axis(TR*basereflection)
        reflectPose_f = - cross(vec(axis_f), vec(axis_t))
        if dot(reflectPose_f, reflectPose) < 0
            reflectPose_f = -reflectPose_f
        end
        reflectPose_f = reflectPose_f ./ norm(reflectPose_f)
    else 
        reflectPose_f = reflectPose
    end
    vec(axis_f), vec(reflectPose_f), degree_f, isreflect
end

@everywhere function refineAxis_D(points, axis, degree, reflectPose, log)
    @printf(log, "%s\n", "refine_D:")
    axis, reflectPose, degree, isreflect = refineAxis_C(points, axis, degree, reflectPose, log)
    symType = ""
    if degree == 1 && isreflect
        axis_t = refineAxis_reflect(points, axis, log)
        if isempty(axis_t) 
            symType = "Cs"
        else
            symType = "C2v"
            degree = 2
            axis = reflectPose
            reflectPose = axis_t
        end
        return vec(axis), vec(reflectPose), degree, symType
    elseif degree == 1
        axis_t = refineAxis_reflect(points, axis, log)
        if isempty(axis_t)
            symType = "E"
            axis_f = axis
            reflectPose_f = reflectPose
        else
            symType = "Cs"
            axis_f = cross(vec(axis_t), vec(reflectPose))
            axis_f = axis_f ./ norm(axis_f)
            reflectPose_f = cross(axis_t, axis_f)
        end
        return vec(axis_f), vec(reflectPose_f), degree, symType
    end
    @printf(log, "%s\n", "H plane varify:")
    transform = axis2matrix(axis,-1)
    points_trans = transform * points'
    TR,TT,ER,t = icp(points', points_trans, 1, Minimize="point")
    @printf(log, "%f\n", ER[1])
    if ER[1] < thresh_end()
        if isreflect
            symType = "D$(degree)h"
        else
            symType = "C$(degree)h"
        end
    else
        @printf(log, "%s\n", "degree 2 axises verify:")
        transform = axis2matrix(axis, degree)
        halftransform = axis2matrix(axis, 2*degree)
        baseDRotate = axis2matrix(reflectPose,2)
        if isreflect
            baseDRotate = halftransform * baseDRotate
        end
        valid_D = zeros(degree,1)
        for j = 1:degree
            points_trans = baseDRotate * points'
            TR,TT,ER,t = icp(points', points_trans, 1, Minimize="point")
            @printf(log, "%f ", ER[1])
            if ER[1] < thresh_start()
                TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point")
                @printf(log, "%f", ER[end])
                if ER[end] < thresh_end() valid_D[j]=1; end
            end
            @printf(log, "\n")
            baseDRotate = transform * baseDRotate
        end
        if sum(valid_D)==degree
            if isreflect
                symType = "D$(degree)d"
            else
                symType = "D$(degree)"
            end
        elseif isreflect
            symType = "C$(degree)v"
        else
            symType = "C$(degree)"
        end
    end

    vec(axis), vec(reflectPose), degree, symType
end

@everywhere function estimateDegree(points, axis)
    rbin = 6
    hbin = 4
    abin = 64
    Descriptor = zeros(1, rbin*2*hbin)
    ez = vec(axis) ./ norm(axis)
    ex = cross(Float64[0,0,1], ez)
    if norm(ex)<0.001
        ex = Float64[1,0,0]
    else
        ex = ex ./ norm(ex)
    end
    ey = cross(ez, ex)
    height = points * ez
    inplane = points - height * ez'
    radi = sqrt(sum(inplane.^2, 2))
    inplane = broadcast(./, inplane, radi)
    ang = acos(clamp(inplane * ex, -1.0, 1.0))  #angle is another function
    reflect = inplane * ey
    f_reflect(x) = x<0
    idx = find(f_reflect, reflect)
    ang[idx] = broadcast(+, 2*pi, -ang[idx])
    height_q = floor(Int32, height .* hbin) + 1 + hbin
    height_q = min(height_q, 2*hbin)
    height_q = max(height_q, 1)
    radi_q = floor(Int32, radi .* rbin) + 1
    radi_q = min(radi_q, rbin)
    angle_q = floor(Int32, ang ./ 2 ./ pi .* abin) + 1
    angle_q = min(angle_q, abin)
    #hist = zeros(abin, 2*hbin, rbin)
    rmax = zeros(abin, 2*hbin)
    for i in 1:size(height_q, 1)
    #    hist[angle_q[i], height_q[i], radi_q[i]]+=1
        rmax[angle_q[i], height_q[i]] = max(radi[i], rmax[angle_q[i], height_q[i]])
    end
    #hist = hist ./ size(points, 1)
    #Fhist = fft(hist, 1)
    #hist = abs(Fhist)
    Frmax = fft(rmax, 1)
    rmax = abs(Frmax)
    rphi = zeros(abin, 2*hbin)
    rphi = angle(Frmax)
    maxval, Descriptor = findmax(rmax[2:end,:],1)
    #maxval, Descriptor = findmax(hist[2:end,:,:],1)
    Descriptor = (Descriptor-1) .% (abin-1) +1
    Descriptor[maxval .< 0.55] = 20
    phase = [rphi[Descriptor[i]+1, i] for i in 1:2*hbin]
    phase = - vec(phase) ./ vec(Descriptor) + pi/abin
    maxpos = zeros(3, size(phase,1))
    for i = 1:size(phase,1)
        maxpos[:,i] = cos(phase[i]) * ex + sin(phase[i]) * ey
    end
    Descriptor, maxpos
    #hist
end

@everywhere function D2Descriptor(points, numPair, lb=0, ub=2, numBin=40) 
    distance=[]
    numPoints = size(points,1)

    for i = 1:numPair
        v1 = rand(1:numPoints)
        v2 = rand(1:numPoints)
        d = norm(points[v2,:] - points[v1,:])
        push!(distance,d)
    end

    interval = (ub - lb) / numBin + 1e-6
    histogram = zeros(numBin,1)
    for d in distance
        idx = clamp(floor(Int, d/interval)+1, 1, numBin)
        histogram[idx]+=1
    end
    histogram /= numPair
    histogram
end

@everywhere function symLevel(symType)
    if symType=="E"
        return 1
    elseif symType=="Cs"
        return 2
    elseif symType[1]=='C'
        if symType[end]=='h' || symType[end]=='v'
            return 4
        else return 3
        end
    elseif symType[1]=='D'
        if symType[end]=='h'
            return 7
        elseif symType[end]=='d'
            return 6
        else return 5
        end
    else return 8 
    end
end 

@everywhere function detectSelfSymmetry(Mesh, log)

    #mesh PCA and normalization
    face = Mesh.faces;
    Cov = zeros(3,3)
    mcenter = zeros(3,1)
    totalarea = 0
    for f in face
        pt1 = vec(Mesh.vertices[f[1],:])
        pt2 = vec(Mesh.vertices[f[2],:])
        pt3 = vec(Mesh.vertices[f[3],:])
        area = norm(cross(pt1-pt2,pt2-pt3))/2
        Cov = Cov + (pt1*pt1' + pt2*pt2' + pt3*pt3' + 0.5*(pt1*pt2' + pt1*pt3' + pt2*pt1' + pt2*pt3' + pt3*pt1' + pt3*pt2')).*area./6
        mcenter = mcenter + area / 3 * (pt1+pt2+pt3)
        totalarea = totalarea + area
    end
    mcenter = mcenter ./ totalarea
    Cov = (Cov ./ totalarea) - mcenter * mcenter'
    eigenval, eigenvec = eig(Cov)
    idx = sortperm(eigenval)
    eigenval = eigenval[idx]
    eigenvec = eigenvec[:,idx]
    Mesh.vertices = Mesh.vertices .- mcenter'
    diag = maximum(sum(.^(Mesh.vertices,2),2),1)
    Mesh.vertices = Mesh.vertices ./ sqrt(diag)
    eigenval = eigenval ./ (diag) 
    @printf(log, "%s\n", "Mesh PCA-eigenval-eigenvec(colume)")
    @printf(log, "%f %f %f\n", eigenval...)
    @printf(log, "%f %f %f\n", eigenvec[1,:]...)
    @printf(log, "%f %f %f\n", eigenvec[2,:]...)
    @printf(log, "%f %f %f\n", eigenvec[3,:]...)
    
    #fixed point and PCA
    densityProposal = GetSampleDensityProposal(Mesh, 8000)
    densepoints = SamplePoints(Mesh, densityProposal)
    #points = SamplePoints_Voxelize(densepoints, [-1 -1 -1], [1 1 1], 0.025)
    points = densepoints
    desc = @time(ShapeContext(points, 1.5, 6,6,128))
    dists = pairwise(Euclidean(), desc')
    S = exp(- dists.^2 / 0.02)
    Cs = broadcast(./, S, sum(S,2))
    Xs = deepcopy(points)
    for t = 1:10
        Xs = Cs * Xs
    end
    #pygui(true)
    #figure()
    #scatter3D(points[:,1], points[:,2], points[:,3], s=40, c="b")
    #scatter3D(Xs[:,1], Xs[:,2], Xs[:,3], s=70, c="r")
    Xs_c = mean(Xs, 1)
    Cov = Xs'*Xs./size(Xs, 1) - Xs_c'*Xs_c
    val, dir = eig(Cov)
    idx = sortperm(val)
    val = val[idx]
    dir = dir[:, idx]
    @printf(log, "%s\n", "Fixed point PCA-eigenval-eigenvec(column)")
    @printf(log, "%f %f %f\n", val...)
    @printf(log, "%f %f %f\n", dir[1,:]...)
    @printf(log, "%f %f %f\n", dir[2,:]...)
    @printf(log, "%f %f %f\n", dir[3,:]...)

    #suggest promising symmetry
    symType = "E"
    canonical_dir = eigenvec
    if val[1] > 0.001
        #no symmetry, return
        @printf(log, "%s\n", "no symmetry")
        canonical_dir = eigenvec
        symType = "E"

    elseif val[2] > 0.001
        #reflection normal = first eigenvector
        axis = dir[:,1]
        degree = -1
        axis_f = refineAxis_reflect(densepoints, axis, log)
        if isempty(axis_f)
            @printf(log, "%s\n", "reflection test failed")
            canonical_dir = eigenvec
            symType = "E"

        else
            @printf(log, "%s\n", "reflection test pass")
            axis2 = cross(dir[:,3], axis_f)
            axis2 = axis2 ./ norm(axis2)
            axis3 = cross(axis_f, axis2)
            Cov_2 = vcat(axis2', axis3') * Cov * hcat(axis2, axis3)
            val2,dir2 = eig(Cov_2)
            dir2 = dir2[:, sortperm(val2)]
            axis2_f=hcat(axis2, axis3)*dir2[:,1] 
            axis3_f=hcat(axis2, axis3)*dir2[:,2]
            canonical_dir = hcat(axis_f, axis2_f, axis3_f)
            symType = "Cs"
        end
    elseif val[3] > 0.001
        #rotation axis = third eigenvector
        axis = dir[:,3]
        @printf(log, "%s\n", "C group axis")
        degrees, reflectPoses = estimateDegree(densepoints, axis)
        @printf(log, "%f %f %f %f %f %f %f %f\n", degrees...)
        idx = sortperm(vec(degrees), rev=true)
        degrees = degrees[idx]
        reflectPoses = reflectPoses[:,idx]
        degree_f = 0
        axis_f = []
        reflectPose = []
        isreflect_f = false
        for i = 1:size(degrees,1)
            if (i>1) && (degrees[i]==degrees[i-1]) continue; end
            axis_t, reflectPose_t, degree_t, isreflect= refineAxis_C(densepoints, axis, degrees[i], reflectPoses[:,i], log)
            if degree_t > degree_f 
                degree_f = degree_t
                axis_f = axis_t
                reflectPose = reflectPose_t
                isreflect_f = isreflect
            elseif degree_t == degree_f && isreflect
                axis_f = axis_t
                reflectPose = reflectPose_t
                isreflect_f = isreflect
            end
        end
        canonical_dir[:,2] = axis_f
        canonical_dir[:,3] = reflectPose
        canonical_dir[:,3] -= canonical_dir[:,2]*dot(canonical_dir[:,2],canonical_dir[:,3])
        canonical_dir[:,3] /= norm(canonical_dir[:,3])
        canonical_dir[:,1] = cross(canonical_dir[:,2],canonical_dir[:,3]) 
        symType = "C$(degree_f)"
        if isreflect_f
            symType = symType * "v"
        end
        if degree_f==1
            symType = "Cs"
            if !isreflect_f
                symType = "E"
            end
        end
    elseif val[3] <= 0.001
        # higher order symmetry
        if eigenval[1]/eigenval[2] < 0.9 && eigenval[2]/eigenval[3] < 0.9
            @printf(log, "%s\n", "cuboid group")
            degrees = ones(3,1)
            reflect = falses(3,1)
            reflectPoses = zeros(3,3)
            for a = 1:3
                axis = eigenvec[:,a]
                if a<3
                    reflectPose = eigenvec[:,a+1]
                else
                    reflectPose = eigenvec[:,1]
                end
                canonical_dir[:,a], reflectPose_t, degree_t, isreflect_t = refineAxis_C(densepoints, axis, 2, reflectPose, log)
                if degree_t == 2 degrees[a]=2; end
                if isreflect_t reflect[a]=1; reflectPoses[:,a] = reflectPose_t; end
            end
            @printf(log, "%d %d %d\n", degrees...)
            @printf(log, "%d %d %d\n", reflect...)
            idx = find(degrees.==2)
            if isempty(idx)
                if reflect[1] || reflect[2] || reflect[3]
                    idx_r = find(reflect.==1)
                    if size(idx_r,1)>2
                        @printf(log, "%s\n", "error: two reflection planes must result in C2")
                    else
                        idx_r = idx_r[1]
                        canonical_dir[:,2] = canonical_dir[:,idx_r]
                        canonical_dir[:,3] = reflectPoses[:,idx_r]
                        canonical_dir[:,3] -= canonical_dir[:,2]*dot(canonical_dir[:,2],canonical_dir[:,3])
                        canonical_dir[:,3] /= norm(canonical_dir[:,3])
                        canonical_dir[:,1] = cross(canonical_dir[:,2],canonical_dir[:,3])
                        symType = "Cs"
                    end
                else 
                    canonical_dir = eigenvec
                    symType = "E"
                end    
            elseif size(idx,1)==3
                if reflect[1] && reflect[2] && reflect[3]
                    symType = "D2h"
                elseif !(reflect[1] || reflect[2] || reflect[3])
                    symType = "D2"
                else 
                    @printf(log, "%s\n", "error: cuboid groups should not have D2d")
                end
                canonical_dir[:,2]=canonical_dir[:,2] - canonical_dir[:,3]*dot(canonical_dir[:,2],canonical_dir[:,3])
                canonical_dir[:,2]=canonical_dir[:,2] ./ norm(canonical_dir[:,2])
                canonical_dir[:,1]=cross(canonical_dir[:,2], canonical_dir[:,3])
            elseif size(idx,1)==1
                idx = idx[1]
                if reflect[idx]
                    symType = "C2v"
                else
                    axis_t = refineAxis_reflect(densepoints, canonical_dir[:,idx], log)
                    if isempty(axis_t)
                        symType = "C2"
                    else
                        symType = "C2h"
                        canonical_dir[:,idx]=axis_t
                    end
                end
                temp = canonical_dir[:,idx]
                canonical_dir[:,idx]=canonical_dir[:,2]
                canonical_dir[:,2]=temp
                canonical_dir[:,3]=canonical_dir[:,3] - canonical_dir[:,2]*dot(canonical_dir[:,2],canonical_dir[:,3])
                canonical_dir[:,3]=canonical_dir[:,3] ./ norm(canonical_dir[:,3])
                canonical_dir[:,1]=cross(canonical_dir[:,2], canonical_dir[:,3])
            else 
                @printf(log, "%s\n", "error resolving symType C2/D2")
            end
        elseif eigenval[1]/eigenval[2] < 0.9 || eigenval[2]/eigenval[3] < 0.9
            @printf(log, "%s\n", "D group axis")
            if eigenval[1]/eigenval[2] < 0.9
                axis = eigenvec[:,1]
            else 
                axis = eigenvec[:,3]
            end
            degrees, reflectPoses = estimateDegree(densepoints, axis)
            @printf(log, "%f %f %f %f %f %f %f %f\n", degrees...)
            idx = sortperm(vec(degrees), rev=true)
            degrees = degrees[idx]
            reflectPoses = reflectPoses[:,idx]
            degree_f = 0
            axis_f = []
            reflectPose = []
            symType_t = ""
            for i = 1:size(degrees,1)
                if (i>1) && (degrees[i]==degrees[i-1]) continue; end
                axis_t, reflectPose_t, degree_t, symType_t= refineAxis_D(densepoints, axis, degrees[i], reflectPoses[:,i], log)
                if degree_t > degree_f 
                    degree_f = degree_t
                    axis_f = axis_t
                    reflectPose = reflectPose_t
                    symType = symType_t
                elseif degree_t == degree_f && symLevel(symType_t) > symLevel(symType)
                    axis_f = axis_t
                    reflectPose = reflectPose_t
                    symType = symType_t
                end
            end
            canonical_dir[:,2] = axis_f
            canonical_dir[:,3] = reflectPose
            canonical_dir[:,3] -= canonical_dir[:,2]*dot(canonical_dir[:,2],canonical_dir[:,3])
            canonical_dir[:,3] /= norm(canonical_dir[:,3])
            canonical_dir[:,1] = cross(canonical_dir[:,2],canonical_dir[:,3]) 

        else
            @printf(log, "%s\n", "calling symSpectral")
            candtrans = symmetrySpectral_short(points,desc)
            symtrans, symscores = refineTransform(candtrans, densepoints)
            axises_f, degrees_f, reflections_f = processSymmetry(symtrans, symscores)
            num_high_degree = 0
            highest_deg = 1
            if isempty(axises_f)
                @printf(log, "%s\n", "no symmetry found")
            else
                @printf(log, "%s\n", "axises found:")
                for i = 1:size(axises_f, 2)
                    @printf(log, "%f %f %f\n", axises_f[:,i]...)  
                    @printf(log, "degree %d\n", degrees_f[i])
                    @printf(log, "reflect %d\n", reflections_f[i])
                end
                for i = 1:size(axises_f, 2)
                    if (degrees_f[i] != 1)
                        degrees, reflectPoses = estimateDegree(densepoints, axises_f[:,i])
                        idx = sortperm(vec(degrees), rev=true)
                        degrees = degrees[idx]
                        reflectPoses = reflectPoses[:,idx]  
                        degree_f = 0
                        symType_f = "E"
                        axis_f = []
                        reflectPose = []
                        symType_t = ""
                        for idx = 1:size(degrees,1)
                            if (idx>1) && (degrees[idx]==degrees[idx-1]) continue; end
                            axis_t, reflectPose_t, degree_t, symType_t= refineAxis_D(densepoints, axises_f[:,i], degrees[idx], reflectPoses[:,idx], log)
                            if degree_t > degree_f
                                degree_f = degree_t
                                axis_f = axis_t
                                reflectPose = reflectPose_t
                                symType_f = symType_t
                            elseif degree_t == degree_f && symLevel(symType_t) > symLevel(symType_f)
                                axis_f = axis_t
                                reflectPose = reflectPose_t
                                symType_f = symType_t
                            end
                        end
                        if degree_f > 2
                            num_high_degree += 1
                        end
                        println(degree_f)
                        println(axis_f)
                        if (degree_f > highest_deg) || (degree_f == highest_deg && symLevel(symType_f) > symLevel(symType)) 
                            highest_deg = degree_f
                            symType = symType_f
                            canonical_dir[:,2] = axis_f
                            canonical_dir[:,3] = reflectPose
                            canonical_dir[:,3] -= canonical_dir[:,2]*dot(canonical_dir[:,2],canonical_dir[:,3])
                            canonical_dir[:,3] /= norm(canonical_dir[:,3])
                            canonical_dir[:,1] = cross(canonical_dir[:,2],canonical_dir[:,3])
                            println(canonical_dir)
                        elseif (degree_f == highest_deg) && (symLevel(symType_f) == symLevel(symType))
                            if (highest_deg == 5) && dot(axis_f, canonical_dir[:,2]) < 0
                                axis_f = -axis_f
                            end
                            canonical_dir[:,3] = axis_f - canonical_dir[:,2] .* dot(canonical_dir[:,2], axis_f)
                            canonical_dir[:,3] /= norm(canonical_dir[:,3])
                            canonical_dir[:,1] = cross(canonical_dir[:,2],canonical_dir[:,3])

                        end
                    end
                end
                if (num_high_degree > 1)
                    println("special group")
                    if highest_deg == 5
                        symType = "I"
                    elseif highest_deg == 4
                        symType = "O"
                    elseif highest_deg == 3 
                        symType = "T"
                    else
                        symType = "O3"
                        canonical_dir = eye(3)
                    end
                end
                println(symType)
                println(canonical_dir)
            end                
        end
    end
    points_new = points * canonical_dir
    thirdMoment = sum(points_new.^3, 1)
    for m = 1:3
        if thirdMoment[m]<0 && symLevel(symType) < 6
            canonical_dir[:,m]=-canonical_dir[:,m]
        end
    end

    symType, canonical_dir, mcenter
end


@everywhere function main(synsetID, modelname)
    println(modelname)
    model = load("/orions3-zfs/projects/haosu/ShapeNetCore2015Spring/ShapeNetCore.v1/" * synsetID * "/" * modelname * "/model.obj")
    v= zeros(size(model.vertices,1),3)
    for i = 1:size(model.vertices,1)
        v[i,:]=[model.vertices[i][1], model.vertices[i][2], model.vertices[i][3]]
    end
    f = zeros(size(model.faces, 1),3)
    for i = 1:size(model.faces,1)
        f[i,:] = [model.faces[i][1]+1, model.faces[i][2]+1, model.faces[i][3]+1]
    end
    newMesh = BuildMesh(v,f)
    vmin = minimum(newMesh.vertices, 1)
    vmax = maximum(newMesh.vertices, 1)

    center = ((vmax+vmin) * 0.5)
    radius = norm(vmax-vmin) * 0.5
    newMesh.vertices = (newMesh.vertices .- center) ./ radius;
    logname = "Results/" * synsetID * "/" * modelname * ".log"
    fout = open(logname, "w")
    symType, canonical, mcenter = detectSelfSymmetry(newMesh, fout)
    close(fout)
    filename = "Results/" * synsetID * "/" * modelname * ".sym3"
    fout = open(filename, "w")
    @printf(fout, "%s\n", symType)
    @printf(fout, "%f %f %f\n", canonical[1,:]...)
    @printf(fout, "%f %f %f\n", canonical[2,:]...)
    @printf(fout, "%f %f %f\n", canonical[3,:]...)
    close(fout)
end



#synsetID = "02942699" 
#println(synsetID)
#models = readall("/orions4-zfs/projects/haohe/3DSIChallenge/deduplicate/deduplicated_model_list/" * synsetID * ".txt")
models = readall("tasklist.txt")
models = split(models, '\n')
models = [split(models[i], '/') for i = 1:size(models,1)-1 ]
synsets = [models[i][1] for i = 1:size(models,1)]
models = [models[i][2] for i = 1:size(models,1)]
#for i = size(models,1):-1:1
#    if isfile("Results/"* synsets[i] * "/" * models[i] * ".sym3") && isfile("Results/"*synsets[i]*"/"*models[i]*".log")
#        deleteat!(models, i)
#        deleteat!(synsets, i)
#    end
#end
#println(size(models))
#shuffle!(models)
#synsets = fill(synsetID, size(models))
#for i = 1:size(models,1)
#    main(synsetID, models[i])
#end
pmap(main, synsets, models)


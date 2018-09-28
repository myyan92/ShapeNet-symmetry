push!(LOAD_PATH, pwd())

# numCore = 10
# addprocs(numCore - 1)

@everywhere using MeshIO
@everywhere using FileIO
@everywhere using SamplePointsUtil
include("./io.jl")
include("./symmetrySpectral.jl")
include("./icp.jl")
include("./refineAxis.jl")

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
    rmax = zeros(abin, 2*hbin)
    for i in 1:size(height_q, 1)
        rmax[angle_q[i], height_q[i]] = max(radi[i], rmax[angle_q[i], height_q[i]])
    end
    Frmax = fft(rmax, 1)
    rmax = abs(Frmax)
    rphi = zeros(abin, 2*hbin)
    rphi = angle(Frmax)
    maxval, Descriptor = findmax(rmax[2:end,:],1)
    Descriptor = (Descriptor-1) .% (abin-1) +1
    Descriptor[maxval .< 0.55] = 20
    Descriptor[Descriptor .> 20] = 20
    phase = [rphi[Descriptor[i]+1, i] for i in 1:2*hbin]
    phase = - vec(phase) ./ vec(Descriptor) + pi/abin
    maxpos = zeros(3, size(phase,1))
    for i = 1:size(phase,1)
        maxpos[:,i] = cos(phase[i]) * ex + sin(phase[i]) * ey
    end
    Descriptor, maxpos
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

@everywhere function get_canonical(reflectNormal, rotationAxis, reflectPose, symType)
    canonical_dir = eye(3)
    if symType=="E"
        return eye(3)
    elseif symType=="Cs"
        if reflectNormal==nothing
            reflectNormal = cross(rotationAxis, reflectPose)
            reflectNormal /= norm(reflectNormal)
        end
        canonical_dir[:,3] = reflectNormal
        temp = cross(canonical_dir[:,2], canonical_dir[:,3])
        if norm(temp) > 0.2
            canonical_dir[:,1] = temp / norm(temp)
            canonical_dir[:,2] = cross(canonical_dir[:,3], canonical_dir[:,1])
        else
            canonical_dir[:,2] = cross(canonical_dir[:,3], canonical_dir[:,1])
            canonical_dir[:,2] /= norm(canonical_dir[:,2])
            canonical_dir[:,1] = cross(canonical_dir[:,2], canonical_dir[:,3])
        end
    else
        canonical_dir[:,2] = rotationAxis
        canonical_dir[:,1] = reflectPose
        canonical_dir[:,1] -= canonical_dir[:,2]*dot(canonical_dir[:,2],canonical_dir[:,1])
        canonical_dir[:,1] /= norm(canonical_dir[:,1])
        canonical_dir[:,3] = cross(canonical_dir[:,1],canonical_dir[:,2])
    end
    return canonical_dir
end

@everywhere function detectSelfSymmetry(Mesh, log)

    #mesh PCA
    face = Mesh.faces;
    Cov = zeros(3,3)
    totalarea = 0
    for f in face
        pt1 = vec(Mesh.vertices[f[1],:])
        pt2 = vec(Mesh.vertices[f[2],:])
        pt3 = vec(Mesh.vertices[f[3],:])
        area = norm(cross(pt1-pt2,pt2-pt3))/2
        Cov = Cov + (pt1*pt1' + pt2*pt2' + pt3*pt3' + 0.5*(pt1*pt2' + pt1*pt3' + pt2*pt1' + pt2*pt3' + pt3*pt1' + pt3*pt2')).*area./6
        totalarea = totalarea + area
    end
    Cov = Cov ./ totalarea
    eigenval, eigenvec = eig(Cov)
    idx = sortperm(eigenval)
    eigenval = eigenval[idx]
    eigenvec = eigenvec[:,idx]

    @printf(log, "%s\n", "Mesh PCA-eigenval-eigenvec(colume)")
    @printf(log, "%f %f %f\n", eigenval...)
    @printf(log, "%f %f %f\n", eigenvec[1,:]...)
    @printf(log, "%f %f %f\n", eigenvec[2,:]...)
    @printf(log, "%f %f %f\n", eigenvec[3,:]...)
   
    #fixed point and PCA
    densepoints = SamplePoints(Mesh.vertices, Mesh.faces, 8000)
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
    canonical_dir = eye(3)
    translate = zeros(3,1)
    if val[1] > 0.001
        #no symmetry, return
        @printf(log, "%s\n", "no symmetry")

    elseif val[2] > 0.001
        #reflection normal = first eigenvector
        axis = dir[:,1]
        degree = -1
        axis_f, translate_f = refineAxis_reflect(densepoints, axis, log)
        if isempty(axis_f)
            @printf(log, "%s\n", "reflection test failed")
        else
            @printf(log, "%s\n", "reflection test pass")
            translate = translate_f
            symType = "Cs"
        end
        canonical_dir = get_canonical(axis_f, nothing, nothing, symType)
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
        translate_f = []
        isreflect_f = false
        for i = 1:size(degrees,1)
            if (i>1) && (degrees[i]==degrees[i-1]) continue; end
            axis_t, reflectPose_t, translate_t, degree_t, isreflect= refineAxis_C(densepoints, axis, degrees[i], reflectPoses[:,i], log)
            if degree_t > degree_f 
                degree_f = degree_t
                axis_f = axis_t
                reflectPose = reflectPose_t
                translate_f = translate_t
                isreflect_f = isreflect
            elseif degree_t == degree_f && isreflect
                axis_f = axis_t
                reflectPose = reflectPose_t
                translate_f = translate_t
                isreflect_f = isreflect
            end
        end
        translate = translate_f
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
        canonical_dir = get_canonical(nothing, axis_f, reflectPose, symType)
    elseif val[3] <= 0.001
        # higher order symmetry
        if eigenval[1]/eigenval[2] < 0.9 && eigenval[2]/eigenval[3] < 0.9
            @printf(log, "%s\n", "cuboid group")
            axis = eigenvec[:,1]
            reflectPose = eigenvec[:,2]
            axis_f, reflectPose, translate, degree_f, symType = refineAxis_D(densepoints, axis, 2, reflectPose, log)
            canonical_dir = get_canonical(nothing, axis_f, reflectPose, symType)

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
                axis_t, reflectPose_t, translate_t, degree_t, symType_t= refineAxis_D(densepoints, axis, degrees[i], reflectPoses[:,i], log)
                if degree_t > degree_f 
                    degree_f = degree_t
                    axis_f = axis_t
                    reflectPose = reflectPose_t
                    translate = translate_t
                    symType = symType_t
                elseif degree_t == degree_f && symLevel(symType_t) > symLevel(symType)
                    axis_f = axis_t
                    reflectPose = reflectPose_t
                    translate = translate_t
                    symType = symType_t
                end
            end
            canonical_dir = get_canonical(nothing, axis_f, reflectPose, symType)

        else
            @printf(log, "%s\n", "calling symSpectral")
            candtrans = symmetrySpectral(points,desc)
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
                        if degrees_f[i] == 0 degrees_f[i]=Int(20); end
                        reflectPose = []
                        for j = 1:size(axises_f, 2)
                            if (reflections_f[j]==-1) && norm(cross(vec(axises_f[:,i]), vec(axises_f[:,j])))>0.9
                                reflectPose = cross(vec(axises_f[:,i]),vec(axises_f[:,j]))
                                break
                            end
                        end
                        if isempty(reflectPose)
                            reflectPose = cross(vec(axises_f[:,i]), vec([1,0,0]))
                            if norm(reflectPose) < 0.2
                                reflectPose = cross(vec(axises_f[:,i]), vec([0,1,0]))
                            end
                        end
                        reflectPose = reflectPose ./ norm(reflectPose)
                        axis_f, reflectPose_f, translate_f, degree_f, symType_f= refineAxis_D(densepoints, axises_f[:,i], degrees_f[i], reflectPose, log)
                        if degree_f > 2
                            num_high_degree += 1
                        end
                        translate = translate .+ translate_f
                        densepoints = densepoints .+ translate_f'
                        if (degree_f > highest_deg) || (degree_f == highest_deg && symLevel(symType_f) > symLevel(symType)) 
                            highest_deg = degree_f
                            symType = symType_f
                            canonical_dir = get_canonical(nothing, axis_f, reflectPose_f, symType)
                        elseif (degree_f == highest_deg) && (symLevel(symType_f) == symLevel(symType))
                            if (highest_deg == 5) && dot(axis_f, canonical_dir[:,2]) < 0
                                axis_f = -axis_f
                            end
                            canonical_dir[:,1] = axis_f - canonical_dir[:,2] .* dot(canonical_dir[:,2], axis_f)
                            canonical_dir[:,1] /= norm(canonical_dir[:,1])
                            canonical_dir[:,3] = cross(canonical_dir[:,1],canonical_dir[:,2])

                        end
                    end
                end
                if (highest_deg == 1)
                    for i = 1:size(axises_f,2)
                        if reflections_f[i] == -1
                            axis_f, translate_f = refineAxis_reflect(densepoints, axises_f[i], log)
                            if !isempty(axis_f)
                                translate = translate_f
                                symType = "Cs"
                                canonical_dir = get_canonical(axis_f, nothing, nothing, symType)
                            end
                        end
                    end
                end

                if (num_high_degree > 1)
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
            end                
        end
    end
    symType, canonical_dir, translate
end


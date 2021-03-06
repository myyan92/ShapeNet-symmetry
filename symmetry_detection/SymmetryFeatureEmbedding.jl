push!(LOAD_PATH, pwd())
numCore = 10
addprocs(numCore - 1)

@everywhere using MeshIO
@everywhere using FileIO
@everywhere using Distances
@everywhere using SamplePointsUtil
@everywhere using IOUtil
@everywhere using ICPUtil
@everywhere using ShapeContextLib
@everywhere using RefineAxisLib
@everywhere using SymmetrySpectralLib

@everywhere function estimateDegree(points, axis)
    rbin = 6
    hbin = 4
    abin = 64
    Descriptor = zeros(1, rbin*2*hbin)
    ez = vec(axis) / norm(axis)
    ex = cross(Float64[0,0,1], ez)
    if norm(ex)<0.001
        ex = Float64[1,0,0]
    else
        ex = ex / norm(ex)
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
    height_q = floor(Int32, height * hbin) + 1 + hbin
    height_q = min(height_q, 2*hbin)
    height_q = max(height_q, 1)
    radi_q = floor(Int32, radi * rbin) + 1
    radi_q = min(radi_q, rbin)
    angle_q = floor(Int32, ang / 2 / pi * abin) + 1
    angle_q = min(angle_q, abin)
    rmax = zeros(abin, 2*hbin)
    for i in 1:size(height_q, 1)
        rmax[angle_q[i], height_q[i]] = max(radi[i], rmax[angle_q[i], height_q[i]])
    end
    Frmax = fft(rmax, 1)
    rmax = abs(Frmax)
    filter(x) = x>0.8
    idx = find(filter, rmax[1,:])
    rmax = rmax[1:33, idx]
    rphi = zeros(abin, 2*hbin)
    rphi = angle(Frmax)
    rphi = rphi[1:33, idx]
    maxval, Descriptor = findmax(rmax[2:end,:],1)
    Descriptor = (Descriptor-1) .% 32 +1
    for i = 1:size(Descriptor, 2)
        tmpsort = sort(rmax[2:end, i], rev=true)
        tmpidx = sortperm(rmax[2:end, i], rev=true)
        for idx in tmpidx
            if (rmax[idx+1, i] * 1.5 > maxval[i]) && (idx % Descriptor[i] != 0 || Descriptor[i]==1)
                Descriptor[i] = 20
                break
            end
        end
        if (maxval[i] < rmax[1, i] * 0.01) Descriptor[i] = 20 end
    end
    Descriptor[Descriptor .> 20] = 20
    phase = [rphi[Descriptor[i]+1, i] for i in 1:size(Descriptor,1)]
    phase = - broadcast(./, vec(phase), vec(Descriptor)) + pi/abin
    maxpos = zeros(3, size(phase,1))
    for i = 1:size(phase,1)
        maxpos[:,i] = cos(phase[i]) * ex + sin(phase[i]) * ey
    end
 
    Descriptor = vec(Descriptor)
    idx = sortperm(Descriptor, rev=true)
    deduplicate_descriptor = [Descriptor[idx[1]]]
    deduplicate_maxpos = maxpos[:,idx[1]]
    for i = 2:size(Descriptor, 1)
        if Descriptor[idx[i]] != deduplicate_descriptor[end]
            push!(deduplicate_descriptor, Descriptor[idx[i]])
            deduplicate_maxpos = hcat(deduplicate_maxpos, maxpos[:,idx[i]])
        end
    end
    if (size(deduplicate_descriptor, 1)>1) && (deduplicate_descriptor[1]==20)
        deduplicate_descriptor = deduplicate_descriptor[2:end]
        deduplicate_maxpos = deduplicate_maxpos[:,2:end]
    elseif (size(deduplicate_descriptor, 1)==1) && (deduplicate_descriptor[1]==20)
        push!(deduplicate_descriptor, 24)
        deduplicate_maxpos = hcat(deduplicate_maxpos, ex)
    end
    
    deduplicate_descriptor, deduplicate_maxpos
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
    canonical_dir = eye(3,3)
    if symType=="E"
        return canonical_dir
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
        Cov = Cov + (pt1*pt1' + pt2*pt2' + pt3*pt3' + 0.5*(pt1*pt2' + pt1*pt3' + pt2*pt1' + pt2*pt3' + pt3*pt1' + pt3*pt2')) * area / 6
        totalarea = totalarea + area
    end
    Cov = Cov / totalarea
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
    num_points = min(round(Int, 10000/9.0*totalarea), 25000)
    if num_points == 25000
        @printf(log, "large area, cap number of points to 25000.")
    end
    densepoints = SamplePoints(Mesh.vertices, Mesh.faces, num_points) # average distance 0.03
    points = densepoints

    #suggest promising symmetry
    symType = "E"
    canonical_dir = eye(3, 3)
    translate = zeros(3,1)

    if eigenval[1]/eigenval[2] < 0.84 && eigenval[2]/eigenval[3] < 0.84
        @printf(log, "%s\n", "cuboid group")
        axis = eigenvec[:,1]
        reflectPose = eigenvec[:,2]
        axis_f, reflectPose, translate, degree_f, symType = refineAxis_D(densepoints, axis, 2, reflectPose, log)
        canonical_dir = get_canonical(nothing, axis_f, reflectPose, symType)

    elseif eigenval[1]/eigenval[2] < 0.84 || eigenval[2]/eigenval[3] < 0.84
        @printf(log, "%s\n", "D group axis")     
        if eigenval[1]/eigenval[2] < 0.84
            axis = eigenvec[:,1]
        else 
            axis = eigenvec[:,3]
        end
        degrees, reflectPoses = estimateDegree(densepoints, axis)
        for d in degrees
            @printf(log, "%f ", d)
        end
        @printf(log, "\n")
        idx = sortperm(vec(degrees), rev=true)
        degrees = degrees[idx]
        reflectPoses = reflectPoses[:,idx]
        degree_f = 0
        for i = 1:size(degrees,1)
            if (i>1) && (degrees[i]==degrees[i-1]) continue; end
            axis_t, reflectPose_t, translate_t, degree_t, symType_t= refineAxis_D(densepoints, axis, degrees[i], reflectPoses[:,i], log)
            if (degree_t > degree_f) || (degree_t == degree_f && symLevel(symType_t) > symLevel(symType)) 
                degree_f = degree_t
                axis_f = axis_t
                reflectPose = reflectPose_t
                translate = translate_t
                symType = symType_t
            end
        end
        canonical_dir = get_canonical(nothing, axis_f, reflectPose, symType)

    else
        # Requires information from fixed points.
        desc = @time(ShapeContext(points, 1.5, 6,6,128))
        dists = pairwise(Euclidean(), desc')
        # Converting to similarity matrix. In-place to avoid memory usage.
        for i in eachindex(dists)
            dists[i] = exp(- dists[i]^2 / 0.02)
        end
        broadcast!(./, dists, dists, sum(dists, 2))
        Xs = deepcopy(points)
        for t = 1:10
            Xs = dists * Xs
        end
        Xs_c = mean(Xs, 1)
        Cov = Xs'*Xs / size(Xs, 1) - Xs_c'*Xs_c
        val, dir = eig(Cov)
        idx = sortperm(val)
        val = val[idx]
        dir = dir[:, idx]
        @printf(log, "%s\n", "Fixed point PCA-eigenval-eigenvec(column)")
        @printf(log, "%f %f %f\n", val...)
        @printf(log, "%f %f %f\n", dir[1,:]...)
        @printf(log, "%f %f %f\n", dir[2,:]...)
        @printf(log, "%f %f %f\n", dir[3,:]...)

        if val[1] > 0.002
            #no symmetry, return
            @printf(log, "%s\n", "no symmetry")

        elseif val[2] > 0.002
            #reflection normal = first eigenvector
            axis = dir[:,1]
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
            for d in degrees
                @printf(log, "%f ", d)
            end
            @printf(log, "\n")
            idx = sortperm(vec(degrees), rev=true)
            degrees = degrees[idx]
            reflectPoses = reflectPoses[:,idx]
            degree_f = 0
            isreflect_f = false
            for i = 1:size(degrees,1)
                if (i>1) && (degrees[i]==degrees[i-1]) continue; end
                axis_t, reflectPose_t, translate_t, degree_t, symType_t = refineAxis_C(densepoints, axis, degrees[i], reflectPoses[:,i], log)
                if (degree_t > degree_f) || (degree_t == degree_f && symLevel(symType_t) > symLevel(symType))
                    degree_f = degree_t
                    axis_f = axis_t
                    reflectPose = reflectPose_t
                    translate = translate_t
                    symType = symType_t
                end
            end
            canonical_dir = get_canonical(nothing, axis_f, reflectPose, symType)
            
        elseif val[3] <= 0.001
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
                        reflectPose = reflectPose / norm(reflectPose)
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
                            canonical_dir[:,1] = axis_f - canonical_dir[:,2] * dot(canonical_dir[:,2], axis_f)
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
                        canonical_dir = eye(3, 3)
                    end
                end
            end                
        end
    end
    symType, canonical_dir, translate
end


@everywhere function main(synsetID, modelname)
    println(modelname)
    if !isdir("Results/"*synsetID)
        mkpath("Results/"*synsetID)
    end
    newMesh = loadMesh(synsetID, modelname)
    logname = "Results/" * synsetID * "/" * modelname * ".log"
    fout = open(logname, "w")
    symType, canonical, translate = detectSelfSymmetry(newMesh, fout)
    close(fout)
    filename = "Results/" * synsetID * "/" * modelname * ".sym3"
    saveSymmetry(filename, symType, translate, canonical)
end

@everywhere function main_v21(modelname)
    println(modelname)
    filename = "/orions-zfs/projects/mengyuan/export-meshes-simple/meshes/" * modelname * "/" * modelname * ".obj"
    newMesh = loadMesh_v2(filename)
    logname = "results_v21/" * modelname * ".log"
    fout = open(logname, "w")
    symType, canonical, translate = detectSelfSymmetry(newMesh, fout)
    close(fout)
    filename = "results_v21/" * modelname * ".sym3"
    saveSymmetry(filename, symType, translate, canonical)
end


#synsetID = ARGS[1]
#println(synsetID)
#models = readall("./deduplicate_lists/" * synsetID * ".txt")
#models = split(models, '\n')
#synsets = fill(synsetID, size(models,1))
tasklist = ARGS[1]
lines = readlines(open(tasklist))
models = [split(l, '\n')[1] for l in lines]
println(size(models,1))
pmap(main_v21, models)
#lines = readall(tasklist)
#lines = split(lines, '\n')
#pop!(lines)
#synsets = [split(l, ' ')[1] for l in lines]
#models = [split(l, ' ')[2] for l in lines]
#println(size(models,1))
#pmap(main, synsets, models)
#main(synsets[1], models[1])
#main("02747177", "16521a9446e3de14a6f567d4d1e09ecb") # test symSpectral

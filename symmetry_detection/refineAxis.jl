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
    TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point", WorstReject=0.1)
    transform = TR * transform
    @printf(log, "%f %f\n", ER[1], ER[end])
    if ER[1] < thresh_start() && ER[end] < thresh_end()
        axis, angle, reflect = matrix2axis(transform)
        translate = dot(TT,axis)/norm(axis) * axis * (-0.5)
        return axis, translate
    end
    return [], []
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
    translate_f = zeros(3,1)
    if degree_f == 1
        axis_f = axis
    else
        transform = axis2matrix(axis, degree_f)
        baserotation = transform
        except=0
        for j = 1:degree_f-1
            points_trans = baserotation * points'
            TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point", WorstReject=0.1)
            axis_t, degree_t, reflect_t = matrix2axis(TR*baserotation)
            if (degree_t == 0)
                except += 1
            else
                if dot(vec(axis_f), vec(axis_t)) < 0 axis_t = -axis_t; end
                translate = TT - dot(TT, axis_t)/norm(axis_t) * axis_t
                translate = translate + cross(axis_t, translate)./norm(axis_t)./tan(degree_t/2)            
                translate_f -= translate./2
                axis_f = axis_f + axis_t
            end
            baserotation = transform * baserotation
        end
        axis_f = axis_f ./ norm(axis_f)
        translate_f = translate_f ./ (degree_f - except)
    end
    if isreflect;
        transform = axis2matrix(axis_f, degree)
        basereflection = axis2matrix(cross(vec(axis_f),vec(reflectPose)), -1)
        for j = 1:idx-1
            basereflection = transform * basereflection
        end
        points_trans = basereflection * (points' .+ translate_f)
        TR,TT,ER,t = icp(points' .+ translate_f, points_trans, 50, Minimize="point", WorstReject=0.1)
        axis_t, degree_t, reflect_t = matrix2axis(TR*basereflection)
        reflectPose_f = - cross(vec(axis_f), vec(axis_t))
        if dot(reflectPose_f, reflectPose) < 0
            reflectPose_f = -reflectPose_f
        end
        reflectPose_f = reflectPose_f ./ norm(reflectPose_f)
        if (norm(translate_f)==0)
            translate_f = dot(TT, axis_t)/norm(axis_t) .* axis_t .* (-0.5)
        end
    else 
        reflectPose_f = reflectPose
    end
    vec(axis_f), vec(reflectPose_f), vec(translate_f), degree_f, isreflect
end

@everywhere function refineAxis_D(points, axis, degree, reflectPose, log)
    @printf(log, "%s\n", "refine_D:")
    axis, reflectPose, translate, degree, isreflect = refineAxis_C(points, axis, degree, reflectPose, log)
    symType = ""
    if degree == 1 && isreflect
        axis_t, translate_2 = refineAxis_reflect(points, axis, log)
        if isempty(axis_t) 
            symType = "Cs"
        else
            symType = "C2v"
            degree = 2
            axis = reflectPose
            reflectPose = axis_t
            translate +=translate_2
        end
        return vec(axis), vec(reflectPose), vec(translate), degree, symType
    elseif degree == 1
        axis_t, translate = refineAxis_reflect(points, axis, log)
        if isempty(axis_t)
            symType = "E"
            axis_f = axis
            reflectPose_f = reflectPose
            translate = zeros(3,1)
        else
            symType = "Cs"
            axis_f = cross(vec(axis_t), vec(reflectPose))
            axis_f = axis_f ./ norm(axis_f)
            reflectPose_f = cross(axis_t, axis_f)
        end
        return vec(axis_f), vec(reflectPose_f), vec(translate), degree, symType
    end
    @printf(log, "%s\n", "H plane varify:")
    transform = axis2matrix(axis,-1)
    points_trans = transform * points'
    TR,TT,ER,t = icp(points', points_trans, 10, Minimize="point", WorstReject=0.1)
    @printf(log, "%f\n", ER[1])
    if ER[end] < thresh_end()
        translate -= 0.5 * dot(TT, axis) .* axis
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
            points_trans = baseDRotate * (points' .+ translate)
            TR,TT,ER,t = icp(points' .+ translate, points_trans, 1, Minimize="point")
            @printf(log, "%f ", ER[1])
            if ER[1] < thresh_start()
                TR,TT,ER,t = icp(points' .+ translate, points_trans, 50, Minimize="point")
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

    vec(axis), vec(reflectPose), vec(translate), degree, symType
end

type Axis
    class::AbstractString
    degree::Int64
    coordinate::Array{Float64, 2}
end

@everywhere function refineSym(proposal, coordinate, points, log)
    if (proposal.class == "E")
      translate = zeros(3,1)
      symtype = "E"
    elseif (proposal.class == "R")
      axis, translate = refineAxis_reflect(points, coordinate[:,1], log)
      if !isempty(axis)
        symtype = "Cs"
        coordinate[:,1] = axis ./ norm(axis)
        coordinate[:,2] = coordinate[:,2] - dot(coordinate[:,2],coordinate[:,1]) .* coordinate[:,1]
        coordinate[:,2] = coordinate[:,2] ./ norm(coordinate[:,2])
        coordinate[:,3] = cross(coordinate[:,1], coordinate[:,2])
      else
        translate = zeros(3,1)
        symtype = "E"
      end
    elseif (proposal.class == "C")
      axis, reflectPose, translate, degree, isreflect = refineAxis_C(points, coordinate[:,2], proposal.degree, coordinate[:,3], log)
      coordinate[:,2] = axis ./ norm(axis)
      coordinate[:,3] = reflectPose
      coordinate[:,3] = coordinate[:,3] - dot(coordinate[:,3],coordinate[:,2]) .* coordinate[:,2]
      coordinate[:,3] = coordinate[:,3] ./ norm(coordinate[:,3])
      coordinate[:,1] = cross(coordinate[:,2], coordinate[:,3])
      if degree>1
        symtype = "C$degree"
        if (isreflect)
          symtype = symtype * "v"
        end
      elseif isreflect
          symtype = "Cs"
      else
          symtype = "E"
          translate = zeros(3,1)
      end
    elseif (proposal.class == "D")
      axis, reflectPose, translate, degree, symtype = refineAxis_D(points, coordinate[:,2], proposal.degree, coordinate[:,3], log)
      coordinate[:,2] = axis ./ norm(axis)
      coordinate[:,3] = reflectPose
      coordinate[:,3] = coordinate[:,3] - dot(coordinate[:,3],coordinate[:,2]) .* coordinate[:,2]
      coordinate[:,3] = coordinate[:,3] ./ norm(coordinate[:,3])
      coordinate[:,1] = cross(coordinate[:,2], coordinate[:,3])
    end
    return symtype, translate, coordinate
end

@everywhere function loadPC(synsetID, modelname)
    #mesh normalization
    model = load("/orions3-zfs/projects/haosu/ShapeNetCore2015Spring/ShapeNetCore.v1/" * synsetID * "/" * modelname * "/model.obj");
    v= zeros(size(model.vertices,1),3);
    for i = 1:size(model.vertices,1)
        v[i,:]=[model.vertices[i][1], model.vertices[i][2], model.vertices[i][3]];
    end
    f = zeros(size(model.faces, 1),3);
    for i = 1:size(model.faces,1)
        f[i,:] = [model.faces[i][1]+1, model.faces[i][2]+1, model.faces[i][3]+1];
    end
    newMesh = BuildMesh(v,f);

    face = newMesh.faces;
    mcenter = zeros(3,1);
    totalarea = 0;
    for f in face
        pt1 = vec(newMesh.vertices[f[1],:]);
        pt2 = vec(newMesh.vertices[f[2],:]);
        pt3 = vec(newMesh.vertices[f[3],:]);
        area = norm(cross(pt1-pt2,pt2-pt3))/2;
        mcenter = mcenter + area / 3 * (pt1+pt2+pt3);
        totalarea = totalarea + area;
    end
    mcenter = mcenter ./ totalarea;
    newMesh.vertices = newMesh.vertices .- mcenter';
    diag = maximum(sum(.^(newMesh.vertices,2),2),1);
    newMesh.vertices = newMesh.vertices ./ sqrt(diag);

    #sample point cloud
    densityProposal = GetSampleDensityProposal(newMesh, 8000);
    densepoints = SamplePoints(newMesh, densityProposal);
    return densepoints
end


@everywhere function saveSymmetry(filename, symType, translate, coordinate)
    fout = open(filename, "w")
    @printf(fout, "%s\n", symType)
    @printf(fout, "%f %f %f\n", translate[:]...)
    @printf(fout, "%f %f %f\n", coordinate[1,:]...)
    @printf(fout, "%f %f %f\n", coordinate[2,:]...)
    @printf(fout, "%f %f %f\n", coordinate[3,:]...)
    close(fout)
end


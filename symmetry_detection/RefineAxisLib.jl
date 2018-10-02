module RefineAxisLib

push!(LOAD_PATH, pwd())

using MeshIO
using FileIO
using SamplePointsUtil
using ICPUtil
using TransformUtil
#include("./ConnectedComponent.jl")
#include("./symmetrySpectral.jl")

export refineAxis_reflect, refineAxis_C, refineAxis_D, refineSym

function thresh_start()
    return 0.1
end

function thresh_end()
    return 0.03
end

#axis is the normal of the reflection plane. 
#return refined axis or none if not a reflection plane.
function refineAxis_reflect(points, axis, log)
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
function refineAxis_C(points, axis, degree, reflectPose, log)
    @printf(log, "%s\n", "refine_C:")
    transform = axis2matrix(axis, degree)
    baserotation = transform
    basereflection = axis2matrix(cross(axis,reflectPose), -1)
    valid_rotate = zeros(degree,1)
    valid_reflect = zeros(degree,1)
    err_reflect = zeros(degree,1)
    isreflect = true
    for j = 1:degree
        points_trans = baserotation * points'
        TR,TT,ER,t = icp(points', points_trans, 1, Minimize="point")
        @printf(log, "%f ", ER[1])
        if ER[1] < thresh_start()
            TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point")
            @printf(log, "%f", ER[end])
            axis_t, degree_t, reflect_t = matrix2axis(TR*baserotation)
            @printf(log, " %f", degree_t)
            if ER[end] < thresh_end() && (abs(degree_t-6.283/degree*j)<0.15 || abs(6.283-degree_t-6.283/degree*j)<0.15) valid_rotate[j]=1; end
        end
        @printf(log, "\n")
        points_trans = basereflection * points'
        TR,TT,ER,t = icp(points', points_trans, 1, Minimize="point")
        @printf(log, "%f ", ER[1])
        err_reflect[j] = ER[1]
        if ER[1] < thresh_start()
            TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point")
            @printf(log, "%f", ER[end])
            err_reflect[j] = ER[end]
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
            err_reflect = reshape(err_reflect, j, numValid)
            numValid2= sum(valid_reflect, 2)
            errSum= sum(err_reflect, 2)
            minerr = 100
            for kk = 1:j
                if numValid2[kk]==numValid && errSum[kk]<minerr
                    idx = kk
                    minerr = errSum[kk]
                end
            end
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
        for j = 1:degree_f-1
            points_trans = baserotation * points'
            TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point", WorstReject=0.1)
            axis_t, degree_t, reflect_t = matrix2axis(TR*baserotation)
            if degree_t != 0
                if dot(vec(axis), vec(axis_t)) < 0 axis_t = -axis_t; end
                translate = TT - dot(TT, axis_t)/norm(axis_t) * axis_t
                translate = translate + cross(axis_t, translate)./norm(axis_t)./tan(degree_t/2)            
                translate_f -= translate./2
                axis_f = axis_f + axis_t
                points = points .- (translate' ./ 2)
            end
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
        TR,TT,ER,t = icp(points', points_trans, 50, Minimize="point", WorstReject=0.1)
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

function refineAxis_D(points, axis, degree, reflectPose, log)
    @printf(log, "%s\n", "refine_D:")
    axis, reflectPose, translate, degree, isreflect = refineAxis_C(points, axis, degree, reflectPose, log)
    @printf(log, "C degree %d, ", degree)
    if isreflect 
        @printf(log, "%s\n", "reflect")
    end
    axis_t, translate_2 = refineAxis_reflect(points, axis, log)
    if isempty(axis_t)
        @printf(log, "%s\n", "D reflect fail")
    else
        @printf(log, "%s\n", "D reflect pass")
    end
    symType = ""
    degree = min(degree, 20)
    if degree == 1 && isreflect
        if isempty(axis_t) 
            symType = "Cs"
        else
            symType = "C2v"
            degree = 2
            axis = reflectPose
            reflectPose = axis_t
            translate += translate_2
        end
        return vec(axis), vec(reflectPose), vec(translate), degree, symType
    elseif degree == 1
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

    if !isempty(axis_t)
        translate += translate_2
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
        points = points .+ translate'
        for j = 1:degree
            points_trans = baseDRotate * points'
            TR,TT,ER,t = icp(points', points_trans, 1, Minimize="point")
            @printf(log, "%f \n", ER[1])
            if ER[1] < thresh_end() valid_D[j]=1; end
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

function refineSym(proposal, coordinate, points, log)
    if (proposal.class == "E")
      translate = zeros(3,1)
      symtype = "E"
      degree = 1
    elseif (proposal.class == "R")
      axis, translate = refineAxis_reflect(points, coordinate[:,1], log)
      if !isempty(axis)
        symtype = "Cs"
        degree = -1
        coordinate[:,1] = axis ./ norm(axis)
        coordinate[:,2] = coordinate[:,2] - dot(coordinate[:,2],coordinate[:,1]) .* coordinate[:,1]
        coordinate[:,2] = coordinate[:,2] ./ norm(coordinate[:,2])
        coordinate[:,3] = cross(coordinate[:,1], coordinate[:,2])
      else
        translate = zeros(3,1)
        symtype = "E"
        degree = 1
      end
    elseif (proposal.class == "C")
      axis, reflectPose, translate, degree, isreflect = refineAxis_C(points, coordinate[:,2], proposal.degree, coordinate[:,3], log)
      coordinate[:,2] = axis ./ norm(axis)
      coordinate[:,3] = reflectPose
      coordinate[:,3] = coordinate[:,3] - dot(coordinate[:,3],coordinate[:,2]) .* coordinate[:,2]
      coordinate[:,3] = coordinate[:,3] ./ norm(coordinate[:,3])
      coordinate[:,1] = cross(coordinate[:,2], coordinate[:,3])
      degree = min(degree, 20)
      if degree>1
        symtype = "C$degree"
        if (isreflect)
          symtype = symtype * "v"
        end
      elseif isreflect
          symtype = "Cs"
          degree = -1
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
    return symtype, degree, translate, coordinate
end


end # module end

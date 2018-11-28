module RefineAxisLib

push!(LOAD_PATH, pwd())

using MeshIO
using FileIO
using SamplePointsUtil
using ICPUtil
using TransformUtil

export refineAxis_reflect, refineAxis_C, refineAxis_D, refineSym

function thresh_start()
    return 0.15
end

function thresh_end()
    return 0.08
end

#axis is the normal of the reflection plane. 
#return refined axis or none if not a reflection plane.
function refineAxis_reflect(points, axis, log)
    @printf(log, "%s\n", "refine_reflect:")
    transform = axis2matrix(axis, -1)
    points_trans = transform * points'
    TR,TT,ER,HD = icp(points', points_trans, 100, Minimize="point", WorstReject=0.1)
    transform = TR * transform
    @printf(log, "%f %f %f\n", ER[1], ER[end], HD)
    if ER[1] < thresh_start() && HD < thresh_end()
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
    err_rotate = zeros(degree,1)
    valid_reflect = zeros(degree,1)
    err_reflect = zeros(degree,1)
    isreflect = true
    for j = 1:degree
        points_trans = baserotation * points'
        TR,TT,ER,HD = icp(points', points_trans, 1, Minimize="point")
        @printf(log, "%f ", ER[1])
        if ER[1] < thresh_start()
            TR,TT,ER,HD = icp(points', points_trans, 100, Minimize="point")
            @printf(log, "%f %f", ER[end], HD)
            axis_t, degree_t, reflect_t = matrix2axis(TR*baserotation)
            @printf(log, " %f", degree_t)
            if HD < thresh_end() && (abs(degree_t-6.283/degree*j)<0.15 || abs(6.283-degree_t-6.283/degree*j)<0.15) valid_rotate[j]=1; end
        end
        err_rotate[j] = ER[end]
        @printf(log, "\n")
        points_trans = basereflection * points'
        TR,TT,ER,HD = icp(points', points_trans, 1, Minimize="point")
        @printf(log, "%f ", ER[1])
        if ER[1] < thresh_start()
            TR,TT,ER,HD = icp(points', points_trans, 100, Minimize="point")
            @printf(log, "%f %f", ER[end], HD)
            if HD < thresh_end() valid_reflect[j]=1; end  
        end
        err_reflect[j] = ER[end]
        @printf(log, "\n")
        baserotation = transform * baserotation
        basereflection = transform * basereflection
    end

    degree_f = 1
    idx = 0
    for j = 1:degree
        if mod(degree,j) != 0 continue; end
        if all(valid_rotate[j:j:end] .>= 0.5)
            degree_f = div(degree, j)
            break
        end
    end
    
    valid_reflect = reshape(valid_reflect, div(degree, degree_f), degree_f)
    err_reflect = reshape(err_reflect, div(degree, degree_f), degree_f)
    errSum = sum(err_reflect, 2)
    minerr = 100
    for k = 1:size(valid_reflect, 1)
        if all(valid_reflect[k,:] .>= 0.5) && errSum[k] < minerr
            idx = k
            minerr = errSum[k]
        end
    end
    if idx == 0  isreflect = false;  end    

    if degree_f == 2 && (!isreflect)
        if err_rotate[div(degree, 2)] > minimum(err_reflect)
            degree_f = 1
            isreflect = true
            tmpval, idx = findmin(err_reflect)
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
            TR,TT,ER,HD = icp(points', points_trans, 100, Minimize="point", WorstReject=0.1)
            axis_t, degree_t, reflect_t = matrix2axis(TR*baserotation)
            if degree_t != 0
                if dot(vec(axis), vec(axis_t)) < 0 axis_t = -axis_t; end
                translate = TT - dot(TT, axis_t)/norm(axis_t) * axis_t
                translate = translate + cross(axis_t, translate) / norm(axis_t) / tan(degree_t/2)            
                translate_f -= translate / 2
                axis_f = axis_f + axis_t
                points = points .- (translate' / 2)
            end
            baserotation = transform * baserotation
        end
        axis_f = axis_f / norm(axis_f)
    end
    if isreflect
        transform = axis2matrix(axis_f, degree)
        basereflection = axis2matrix(cross(vec(axis_f),vec(reflectPose)), -1)
        for j = 1:idx-1
            basereflection = transform * basereflection
        end
        points_trans = basereflection * points'
        TR,TT,ER,HD = icp(points', points_trans, 100, Minimize="point", WorstReject=0.1)
        axis_t, degree_t, reflect_t = matrix2axis(TR*basereflection)
        reflectPose_f = - cross(vec(axis_f), vec(axis_t))
        if dot(reflectPose_f, reflectPose) < 0
            reflectPose_f = -reflectPose_f
        end
        reflectPose_f = reflectPose_f / norm(reflectPose_f)
        if (norm(translate_f)==0)
            translate_f = dot(TT, axis_t)/norm(axis_t) .* axis_t * (-0.5)
        end
    else 
        reflectPose_f = reflectPose
    end

    symType = "C$(degree_f)"
    if isreflect
        symType = symType * "v"
    end
    if degree_f == 1
        if isreflect
            symType = "Cs"
        else
            symType = "E"
        end
    end
    @printf(log, "C degree %d, ", degree_f)
    if isreflect 
        @printf(log, "%s\n", "reflect")
    end
    vec(axis_f), vec(reflectPose_f), vec(translate_f), degree_f, symType
end

function refineAxis_D(points, axis, degree, reflectPose, log)
    @printf(log, "%s\n", "refine_D:")
    axis, reflectPose, translate, degree, symType = refineAxis_C(points, axis, degree, reflectPose, log)

    axis_h, translate_h = refineAxis_reflect(points, axis, log)
    if isempty(axis_h)
        @printf(log, "%s\n", "D reflect fail")
    else
        @printf(log, "%s\n", "D reflect pass")
    end

    degree = min(degree, 20)
    if symType == "Cs" 
        if !isempty(axis_h) 
            symType = "C2v"
            degree = 2
            axis = reflectPose
            reflectPose = axis_h
            translate += translate_h
        end
        return vec(axis), vec(reflectPose), vec(translate), degree, symType
    elseif symType == "E"
        if !isempty(axis_h)
            symType = "Cs"
            axis = cross(vec(axis_h), vec(reflectPose))
            axis = axis / norm(axis)
            reflectPose = cross(axis_h, axis)
        end
        return vec(axis), vec(reflectPose), vec(translate), degree, symType
    end

    if !isempty(axis_h)
        translate += translate_h
        if symType[end] == 'v'
            symType = "D$(degree)h"
        else
            symType = "C$(degree)h"
        end
        return vec(axis), vec(reflectPose), vec(translate), degree, symType        
    else
        @printf(log, "%s\n", "degree 2 axises verify:")
        transform = axis2matrix(axis, degree)
        halftransform = axis2matrix(axis, 2*degree)
        baseDRotate = axis2matrix(reflectPose,2)
        if symType[end] == 'v'
            baseDRotate = halftransform * baseDRotate
        end
        valid_D = zeros(degree,1)
        points = points .+ translate'
        for j = 1:degree
            points_trans = baseDRotate * points'
            TR,TT,ER,HD = icp(points', points_trans, 1, Minimize="point")
            @printf(log, "%f %f\n", ER[1], HD)
            if HD < thresh_end() valid_D[j]=1; end
            baseDRotate = transform * baseDRotate
        end
        if all(valid_D .>= 0.5)
            if symType[end] == 'v'
                symType = "D$(degree)d"
            else
                symType = "D$(degree)"
            end
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
            coordinate[:,1] = axis / norm(axis)
            coordinate[:,2] = coordinate[:,2] - dot(coordinate[:,2],coordinate[:,1]) .* coordinate[:,1]
            coordinate[:,2] = coordinate[:,2] / norm(coordinate[:,2])
            coordinate[:,3] = cross(coordinate[:,1], coordinate[:,2])
        else
            translate = zeros(3,1)
            symtype = "E"
            degree = 1
        end
    elseif (proposal.class == "C")
        axis, reflectPose, translate, degree, symType = refineAxis_C(points, coordinate[:,2], proposal.degree, coordinate[:,3], log)
        coordinate[:,2] = axis / norm(axis)
        coordinate[:,3] = reflectPose
        coordinate[:,3] = coordinate[:,3] - dot(coordinate[:,3],coordinate[:,2]) .* coordinate[:,2]
        coordinate[:,3] = coordinate[:,3] / norm(coordinate[:,3])
        coordinate[:,1] = cross(coordinate[:,2], coordinate[:,3])
        degree = min(degree, 20)
        if symtype == "Cs"
            degree = -1
        end
    elseif (proposal.class == "D")
        axis, reflectPose, translate, degree, symtype = refineAxis_D(points, coordinate[:,2], proposal.degree, coordinate[:,3], log)
        coordinate[:,2] = axis / norm(axis)
        coordinate[:,3] = reflectPose
        coordinate[:,3] = coordinate[:,3] - dot(coordinate[:,3],coordinate[:,2]) .* coordinate[:,2]
        coordinate[:,3] = coordinate[:,3] / norm(coordinate[:,3])
        coordinate[:,1] = cross(coordinate[:,2], coordinate[:,3])
        degree = min(degree, 20)
        if symtype == "Cs"
            degree = -1
        end
    end
    return symtype, degree, translate, coordinate
end


end # module end

module TransformUtil

export axis2matrix, matrix2axis

function axis2matrix(axis, degree)
    matrix = eye(3)
    if degree == -1
        axis = axis ./ norm(axis)
        matrix = matrix - 2 * axis * axis'
    elseif degree == 1
        matrix = eye(3)
    else 
        angle = 2 * pi / degree
        a = cos(angle/2)
        axis = axis ./ norm(axis) .* sin(angle/2)
        matrix[1,1] = matrix[1,1] - 2 * (axis[2]^2 + axis[3]^2)
        matrix[2,2] = matrix[2,2] - 2 * (axis[1]^2 + axis[3]^2)
        matrix[3,3] = matrix[3,3] - 2 * (axis[1]^2 + axis[2]^2)
        matrix[1,2] = 2 * (axis[1] * axis[2] - axis[3] * a)
        matrix[2,1] = 2 * (axis[1] * axis[2] + axis[3] * a)
        matrix[1,3] = 2 * (axis[1] * axis[3] + axis[2] * a)
        matrix[3,1] = 2 * (axis[1] * axis[3] - axis[2] * a)
        matrix[2,3] = 2 * (axis[2] * axis[3] - axis[1] * a)
        matrix[3,2] = 2 * (axis[2] * axis[3] + axis[1] * a)
    end
    
    matrix
end

function matrix2axis(transform)
    D, V = eig(transform)
    reflection = sign(det(transform))
    if sum((transform - reflection*eye(3)).^2) < 0.001
        axis = [0,0,0]
        ang = 0
    else
        if reflection==1
            tr = transform[1,1] + transform[2,2] + transform[3,3]
            if (tr > 0)  
                S = sqrt(tr+1.0) * 2;
                qw = 0.25 * S;
                qx = (transform[3,2]-transform[2,3]) / S;
                qy = (transform[1,3]-transform[3,1]) / S; 
                qz = (transform[2,1]-transform[1,2]) / S; 
            elseif ((transform[1,1] > transform[2,2])&&(transform[1,1] > transform[3,3]))  
                S = sqrt(1.0 + transform[1,1] - transform[2,2] - transform[3,3]) * 2; 
                qw = (transform[3,2]-transform[2,3]) / S;
                qx = 0.25 * S;
                qy = (transform[2,1]+transform[1,2]) / S; 
                qz = (transform[1,3]+transform[3,1]) / S; 
            elseif (transform[2,2] > transform[3,3])  
                S = sqrt(1.0 + transform[2,2] - transform[1,1] - transform[3,3]) * 2; 
                qw = (transform[1,3]-transform[3,1]) / S;
                qx = (transform[2,1]+transform[1,2]) / S; 
                qy = 0.25 * S;
                qz = (transform[3,2]+transform[2,3]) / S; 
            else  
                S = sqrt(1.0 + transform[3,3] - transform[1,1] - transform[2,2]) * 2; 
                qw = (transform[2,1]-transform[1,2]) / S;
                qx = (transform[1,3]+transform[3,1]) / S;
                qy = (transform[3,2]+transform[2,3]) / S;
                qz = 0.25 * S;
            end
            axis = real([qx,qy,qz]) ./ sqrt(1-qw*qw)
            ang = acos(qw) * 2
        else
            idx = indmin(abs(D-reflection))
            axis = real(V[:,idx])
            ang = 0
        end
    end

    axis, ang, reflection
end    


end # module end

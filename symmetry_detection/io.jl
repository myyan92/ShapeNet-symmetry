@everywhere type MeshType
  vertices::Array
  faces::Array

  MeshType() = new([], [])
  MeshType(vertices, faces) = new(vertives, faces)
end

@everywhere function loadMesh(synsetID, modelname)
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
    newMesh = MeshType()
    newMesh.vertices = v
    newMesh.faces = []
    for i = 1:size(f, 1)
        push!(newMesh.faces, Array{Int}(f[i, :]))
    end

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

    return newMesh
end


@everywhere function loadMesh_v2(filename)
    #mesh normalization
    model = load(filename);
    v= zeros(size(model.vertices,1),3);
    for i = 1:size(model.vertices,1)
        v[i,:]=[model.vertices[i][1], model.vertices[i][2], model.vertices[i][3]];
    end
    f = zeros(size(model.faces, 1),3);
    for i = 1:size(model.faces,1)
        f[i,:] = [model.faces[i][1]+1, model.faces[i][2]+1, model.faces[i][3]+1];
    end
    newMesh = MeshType()
    newMesh.vertices = v
    newMesh.faces = []
    for i = 1:size(f, 1)
        push!(newMesh.faces, Array{Int}(f[i, :]))
    end

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

    return newMesh
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

function readSymmetry(synsetID, md5)
    filename = "./Results/" * synsetID * "/" * md5 * ".sym3";
    f = open(filename);
    sym_type = readline(f);
    sym_type = chomp(sym_type);
    translate = readline(f);
    translate = split(chomp(translate))
    translate = [float(translate[1]), float(translate[2]), float(translate[3])]
    coordinate = readdlm(f);
    close(f);
    return sym_type, translate, coordinate
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

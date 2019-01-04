module CleanMeshUtil

export MeshType, CleanMesh, IsTriangleMesh

type MeshType
  vertices::Array
  faces::Array
  MeshType() = new([], [])
  MeshType(vertices, faces) = new(vertices, faces)
end

function WeldVertices(mesh::MeshType, eps::Float64)
    vts = mesh.vertices
    vertDict = Dict()
    vertVertDict = Dict()
    k = 1
    newVtsList = []
    for i = 1:size(vts, 1)        
        v = floor(vts[i, :] / eps) * eps
        if !haskey(vertDict, v)
            vertDict[v] = k
            vertVertDict[i] = k
            push!(newVtsList, v)
            k += 1
        else
            vertVertDict[i] = vertDict[v]
        end
    end
    
    newVts = zeros(length(newVtsList), 3)
    for i = 1:length(newVtsList)
        newVts[i, :] = newVtsList[i][:]
    end
    
    newFcs = []
    fcs = mesh.faces
    for f in mesh.faces        
        push!(newFcs, [vertVertDict[fv] for fv in f])
    end
    
    MeshType(newVts, newFcs)
end

function RemoveDegenerateFaces(mesh::MeshType)
    newfaces = []
    for f in mesh.faces
        if length(unique(f))>=3
            push!(newfaces, f)
        end
    end
    mesh.faces = newfaces
    mesh
end

function RemoveDoubleFaces(mesh::MeshType)
    newfaces = []
    for f in mesh.faces
        push!(newfaces, sort(f))
    end
    mesh.faces = unique(newfaces)
    mesh
end

function IsTriangleMesh(mesh::MeshType)
    for f in mesh.faces
        if length(f)!=3
            return false
        end
    end
    true
end

function RemoveFlyingVertices(mesh::MeshType)
    allverts = vcat(mesh.faces...)
    allverts = unique(allverts)
    if length(allverts) < length(mesh.vertices)
        allverts = sort(allverts)
        vertVertDict = Dict()
        for i in 1:length(allverts)
            vertVertDict[allverts[i]]=i
        end
        newVts = zeros(length(allverts), 3)
        for i = 1:length(allverts)
            newVts[i, :] = mesh.vertices[allverts[i], :]
        end
        newFcs = mesh.faces
        for i in 1:length(mesh.faces)
            newFcs[i][1]=vertVertDict[newFcs[i][1]]
            newFcs[i][2]=vertVertDict[newFcs[i][2]]
            newFcs[i][3]=vertVertDict[newFcs[i][3]]
        end
        mesh = MeshType(newVts, newFcs)
    end
    mesh
end

function RemoveIllTriangles(mesh::MeshType)
    newFcs = []
    for f in mesh.faces
        pt1 = vec(mesh.vertices[f[1],:]);
        pt2 = vec(mesh.vertices[f[2],:]);
        pt3 = vec(mesh.vertices[f[3],:]);
        area = norm(cross(pt1-pt2,pt2-pt3))/2;
        diameter = max(norm(pt1-pt2), norm(pt1-pt3), norm(pt2-pt3))
        if area > 0.0001*diameter*diameter
            push!(newFcs, f)
        end
    end
    mesh.faces = newFcs
    mesh
end

function CleanMesh(mesh::MeshType)
    mesh = RemoveIllTriangles(mesh)
    mesh = RemoveFlyingVertices(mesh)
    mesh = WeldVertices(mesh, 1e-7)
    mesh = RemoveDegenerateFaces(mesh)
    mesh = RemoveDoubleFaces(mesh)    
    mesh
end


end # module end

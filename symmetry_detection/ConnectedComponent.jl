function FindRoot!(C, v)
    if C[v] != v
        C[v] = FindRoot!(C, C[v])
    end
    C[v]
end

function GetConnectedComponents(mesh)
    faces = mesh.faces
    verts = Vertices(mesh)
    numVert = size(verts, 1)
    C = collect(1:numVert)
    for f in faces
        roots = []
        for v = 1:size(f,1)
            push!(roots, FindRoot!(C, f[v]))
        end
        minv = minimum(roots)
        for v = 1:size(f,1)
            C[roots[v]] = minv
        end
    end
    for i = 1:size(C,1)
        C[i] = C[C[i]]
    end
    numComp = unique(C)
    components = []
    for i in numComp
        vertSet = find(C .== i)
        vertVertDict = Dict()
        for (j, k) in enumerate(vertSet)
            vertVertDict[k] = j
        end
        subVts = zeros(length(vertSet), 3)        
        for j = 1:length(vertSet)
            subVts[j, :] = verts[vertSet[j], :]
        end
        vertSet = Set(vertSet)
        subFcsList = []
        for f in faces
            if f[1] in vertSet
                newface = []
                for j = 1:size(f,1)
                    push!(newface, vertVertDict[f[j]])
                end
                push!(subFcsList, newface)
            end
        end        
        push!(components, MeshType(subVts, subFcsList))
    end    
    components
end
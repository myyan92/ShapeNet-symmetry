push!(LOAD_PATH, pwd())

#numCore = 4 #10
#addprocs(numCore - 1)

using MeshIO
using FileIO
using IOUtil
using SymmetryFeatureEmbeddingLib

const path2obj = "/orions4-zfs/projects/jyau/SymmDetProject/testing/Minhyuk-style-OBJ/"
# const path2results = "/orions4-zfs/projects/anastasiad/ShapeNet-symmetry/Results/03001627/"
# const path2obj = "/orions3-zfs/projects/anastasiad/ShapeNetSymm/Minhyuk-style-OBJ/"
const path2results = "/orions3-zfs/projects/anastasiad/ShapeNetSymm/Results/03001627/"

@everywhere function main(modelname)
	for partnum = 4 #0:1000
	    filename = "/orions4-zfs/projects/jyau/SymmDetProject/testing/Minhyuk-style-OBJ/" * modelname * "/" * "$partnum.obj"
	    if (!isfile(filename))
	        break
	    end
	    println(filename)
	    newMesh = loadMesh_v2(filename)
            # if (length(newMesh.vertices)<5)
            #     println("mesh has < 5 vertices, continue")
            #     continue
            # end
	    logname = "/orions3-zfs/projects/anastasiad/ShapeNetSymm/Results/03001627/" * modelname * "_$partnum.log"
	    #println(logname)
	    fout = open(logname, "w")
	    symType, canonical, translate, center = detectSelfSymmetry(newMesh, fout)
	    if (symType == "None")
	        println("symType is None, continue")
	        continue
	    end
	    close(fout)
	    symname = "/orions3-zfs/projects/anastasiad/ShapeNetSymm/Results/03001627/" * modelname * "_$partnum.sym"
	    #println(symname)
	    saveSymmetry(symname, symType, translate, canonical, center)
    end
end

# synsetID = ARGS[1]
synsetID = "03001627"
println(synsetID)
# models = readall("./deduplicate_lists/" * synsetID * ".txt")
models = readall("./all_model_lists/" * synsetID * ".txt")
models = split(models, '\n')
println(size(models,1))
for model_id = 625 #:size(models,1)
        println(model_id)
	main(models[model_id])
end
#pmap(main, models)


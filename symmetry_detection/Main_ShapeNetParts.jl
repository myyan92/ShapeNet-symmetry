push!(LOAD_PATH, pwd())

numCore = 10
addprocs(numCore - 1)

using MeshIO
using FileIO
using IOUtil
using SymmetryFeatureEmbeddingLib

const path2obj = "/orions4-zfs/projects/jyau/SymmDetProject/testing/Minhyuk-style-OBJ/"
# const path2results = "/orions4-zfs/projects/anastasiad/ShapeNet-symmetry/Results/03001627/"
# const path2obj = "/orions3-zfs/projects/anastasiad/ShapeNetSymm/Minhyuk-style-OBJ/"
const path2results = "/orions3-zfs/projects/anastasiad/ShapeNetSymm/Results/03001627/"

@everywhere function main(path2obj, modelname)
	for partnum = 0:1000
		filename = "/orions4-zfs/projects/jyau/SymmDetProject/testing/Minhyuk-style-OBJ/" * modelname * "/" * "$partnum.obj"
        if (!isfile(filename))
        	break
        end
	    println(filename)
	    newMesh = loadMesh_v2(filename)
	    logname = "/orions3-zfs/projects/anastasiad/ShapeNetSymm/Results/03001627/" * modelname * "_$partnum.log"
	    println(logname)
	    fout = open(logname, "w")
	    symType, canonical, translate = detectSelfSymmetry(newMesh, fout)
	    close(fout)
	    symname = "/orions3-zfs/projects/anastasiad/ShapeNetSymm/Results/03001627/" * modelname * "_$partnum.sym"
	    println(symname)
	    saveSymmetry(symname, symType, translate, canonical)
    end
end

# synsetID = ARGS[1]
synsetID = "03001627"
println(synsetID)
# models = readall("./deduplicate_lists/" * synsetID * ".txt")
models = readall("./all_model_lists/" * synsetID * ".txt")
models = split(models, '\n')
println(size(models,1))
for model_id = 1:size(models,1)
	main(path2obj[1],models[model_id])
end

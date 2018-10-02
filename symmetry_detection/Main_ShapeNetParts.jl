push!(LOAD_PATH, pwd())

numCore = 10
addprocs(numCore - 1)

@everywhere using MeshIO
@everywhere using FileIO
include("./io.jl")
include("./symmetryFeatureEmbedding_module.jl")

const path2obj = "/orions4-zfs/projects/jyau/SymmDetProject/testing/Minhyuk-style-OBJ/"
const path2results = "/orions4-zfs/projects/anastasiad/ShapeNet-symmetry/Results/03001627/"

@everywhere function main(path2obj, modelname)
	for partnum = 0:50
		filename = "/orions4-zfs/projects/jyau/SymmDetProject/testing/Minhyuk-style-OBJ/" * modelname * "/" * "$partnum.obj"
        if (!isfile(filename))
        	break
        end
	    println(filename)
	    newMesh = loadMesh_v2(filename)
	    logname = "/orions4-zfs/projects/anastasiad/ShapeNet-symmetry/Results/03001627/" * modelname * "_$partnum.log"
	    println(logname)
	    fout = open(logname, "w")
	    symType, canonical, translate = detectSelfSymmetry(newMesh, fout)
	    close(fout)
	    symname = "/orions4-zfs/projects/anastasiad/ShapeNet-symmetry/Results/03001627/" * modelname * "_$partnum.sym"
	    println(symname)
	    saveSymmetry(symname, symType, translate, canonical)
    end
end

# synsetID = ARGS[1]
synsetID = "03001627"
println(synsetID)
models = readall("./deduplicate_lists/" * synsetID * ".txt")
models = split(models, '\n')
#path2obj = fill(path2obj, size(models,1))
println(size(models,1))
# pmap(main, path2obj, models)
main(path2obj[1],models[1])

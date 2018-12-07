push!(LOAD_PATH, pwd())
# numCore = 10
# addprocs(numCore - 1)

using MeshIO
using FileIO
using IOUtil
using SymmetryFeatureEmbeddingLib

@everywhere function main(modelname)
	objname = "/afs/cs.stanford.edu/u/anastasiad/Lingxiao/repptforsiemensfridaymeeting/" * modelname * ".obj"
    println(objname)
    newMesh = loadMesh_v2(objname)
    logname = "/afs/cs.stanford.edu/u/anastasiad/Lingxiao/repptforsiemensfridaymeeting/Results/" * modelname * ".log"
    fout = open(logname, "w")
    symType, canonical, translate = detectSelfSymmetry(newMesh, fout)
    close(fout)
    filename = "/afs/cs.stanford.edu/u/anastasiad/Lingxiao/repptforsiemensfridaymeeting/Results/" * modelname * ".sym"
    saveSymmetry(filename, symType, translate, canonical)
end

tasklist = ARGS[1]
models = readall(tasklist)
models = split(models, '\n')
println(models)
println(size(models,1))
for model_id = 1:size(models,1)
    println(model_id)
	main(models[model_id])
end

push!(LOAD_PATH, pwd())

numCore = 10
addprocs(numCore - 1)

@everywhere using MeshIO
@everywhere using FileIO
include("./io.jl")
include("./symmetryFeatureEmbedding_module.jl")


@everywhere function main(synsetID, modelname)
    println(modelname)
    newMesh = loadMesh(synsetID, modelname)
    logname = "Results/" * synsetID * "/" * modelname * ".log"
    fout = open(logname, "w")
    symType, canonical, translate = detectSelfSymmetry(newMesh, fout)
    close(fout)
    filename = "Results/" * synsetID * "/" * modelname * ".sym3"
    saveSymmetry(filename, symType, translate, canonical)
end

synsetID = ARGS[1]
println(synsetID)
models = readall("./deduplicate_lists/" * synsetID * ".txt")
models = split(models, '\n')
synsets = fill(synsetID, size(models,1))
println(size(models,1))
# pmap(main, synsets, models)
main(synsets[1],models[1])

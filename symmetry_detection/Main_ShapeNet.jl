push!(LOAD_PATH, pwd())
numCore = 10
addprocs(numCore - 1)


using FileIO
using IOUtil
using SymmetryFeatureEmbeddingLib

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

#synsetID = ARGS[1]
#println(synsetID)
#models = read("./deduplicate_lists/" * synsetID * ".txt")
#models = split(String(models), "\n")
#synsets = fill(synsetID, size(models,1))
tasklist = ARGS[1]
lines = read(tasklist)
lines = split(lines, '\n')
pop!(lines)
synsets = [split(l, ' ')[1] for l in lines]
models = [split(l, ' ')[2] for l in lines]
println(size(models,1))
pmap(main, synsets, models)
#main(synsets[1], models[1])
#main("02747177", "16521a9446e3de14a6f567d4d1e09ecb")

#numCore = 20
#addprocs(numCore - 1)

include("./refineAxis.jl")


@everywhere function readSymmetry(synsetID, md5)
    sym_root = "/orions4-zfs/projects/mengyuan/julia_code/symmetry/Results/" * synsetID * "/";
    f = open(sym_root*md5*".sym3");
    sym_type = readline(f);
    sym_type = chomp(sym_type);
    coordinate = readdlm(f);
    if sym_type=="E"
        return [Axis("E", 1, coordinate)];
    end
    if sym_type=="Cs"
        return [Axis("R", -1, coordinate)];
    end
    if (sym_type[1]=='C' || sym_type[1]=='D')
        degree = tryparse(Int64, sym_type[2:end]);
        if isnull(degree)
            degree = parse(Int64, sym_type[2:end-1]);
        else
            degree = get(degree);
        end
    end
    if sym_type[1]=='C'
        if sym_type[end]=='h'
            return [Axis("D", degree, coordinate)];
        else
            return [Axis("C", degree, coordinate)];
        end
    end
    if sym_type[1]=='D'
        return [Axis("D", degree, coordinate)];
    end
    if sym_type[1]=='O'
        if sym_type=="O3"
            degree = 20
        else
            degree = 4
        end
        axises = [Axis("D", degree, coordinate)]
        coordinate = [coordinate[:,2] coordinate[:,3] coordinate[:,1]]
        push!(axises, Axis("D", degree, coordinate))
        coordinate = [coordinate[:,2] coordinate[:,3] coordinate[:,1]]
        push!(axises, Axis("D", degree, coordinate))
        return axises
    end
end


@everywhere function main(synsetID, modelname)
    println(modelname)
    points = loadPC(synsetID, modelname)
    proposal = readSymmetry(synsetID, modelname)
    logname = "Results/" * synsetID * "/" * modelname * ".log2"
    fout = open(logname, "w")
    if size(proposal,1)==1
        proposal = proposal[1]
        symType, translate, coordinate = refineSym(proposal, proposal.coordinate, points, fout)
    else
        num_high_deg = 0
        highest_deg = 1
        highest_coordinate = proposal[1].coordinate
        highest_symType = "E"
        highest_translate = zeros(3,1)
        for i = 1:size(proposal,1)
            symType, translate, coordinate = refineSym(proposal[i], proposal[i].coordinate, points, fout)
            degree = tryparse(Int64, symType[2:end]);
            if isnull(degree)
                degree = parse(Int64, symType[2:end-1]);
            else
                degree = get(degree);
            end
            if degree > 2
                num_high_deg += 1
            end
            if degree > highest_deg
                highest_deg = degree
                highest_coordinate = coordinate
                highest_symType = symType
            end
            highest_translate += translate
            points = points .+ translate'
        end
        if num_high_deg > 1
            if highest_deg == 4
                symType = "O"
            elseif highest_deg == 5
                symType = "I"
            elseif highest_deg == 3
                symType = "T"
            else
                symType = "O3"
            end
        end
        translate = highest_translate
        coordinate = highest_coordinate
    end
    close(fout)
    filename = "Results/" * synsetID * "/" * modelname * ".sym3t"
    saveSymmetry(filename, symType, translate, coordinate)
end

synsetID = ARGS[1] 
println(synsetID)
models = readall("/orions4-zfs/projects/haohe/3DSIChallenge/deduplicate/deduplicated_model_list/" * synsetID * ".txt")
#models = readall("tasklist.txt")
models = split(models, '\n')
#models = [split(models[i], ' ') for i = 1:size(models,1)-1 ]
#synsets = [models[i][1] for i = 1:size(models,1)]
#models = [models[i][2] for i = 1:size(models,1)]
#for i = size(models,1):-1:1
#    if isfile("Results/"* synsets[i] * "/" * models[i] * ".sym3") && isfile("Results/"*synsets[i]*"/"*models[i]*".log")
#        deleteat!(models, i)
#        deleteat!(synsets, i)
#    end
#end
#println(size(models))
#shuffle!(models)
#synsets = fill(synsetID, size(models))
for i = 1:size(models,1)-1
    main(synsetID, models[i])
end
#pmap(main, synsets, models)


include("refineAxis.jl")

synsetID = ARGS[1];
model_list = "/orions4-zfs/projects/haohe/3DSIChallenge/deduplicate/deduplicated_model_list/"*synsetID*".txt";
models = readall(model_list);
models = split(models, '\n');
size(models,1)

graphFile = "/orions4-zfs/projects/mengyuan/julia_code/symmetry/NNGraph/"*synsetID*".txt";
NNgraph = readdlm(graphFile);
size(NNgraph,1)


function readSymmetry(synsetID, md5)
    sym_root = "/orions4-zfs/projects/mengyuan/julia_code/symmetry/Results/" * synsetID * "/";
    f = open(sym_root*md5*".sym3t");
    sym_type = readline(f);
    sym_type = chomp(sym_type);
    translate = readline(f);
    translate = split(chomp(translate))
    translate = [float(translate[1]), float(translate[2]), float(translate[3])]
    coordinate = readdlm(f);
    if sym_type=="E" 
        return []; 
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
    error(" symmetry not implemented")    
end


function mergeAxis(axisList, axis)
    if isempty(axisList)
        return axis
    end
    if !isempty(axis)
        proposals = deepcopy(axis);
    else
        proposals = [deepcopy(axisList[1])]
    end
    for i = 1:size(axisList,1)
        merged = false
        for j = 1:size(proposals,1)
            if merged continue; end
            merged = true
            vec1 = proposals[j].coordinate[:,2]
            if proposals[j].class=="R"
                vec1 = proposals[j].coordinate[:,1]
            end
            vec2 = axisList[i].coordinate[:,2]
            if axisList[i].class=="R"
                vec2 = axisList[i].coordinate[:,1]
            end
            dist = abs(vec2'*vec1);
            if dist[1] > 0.97
                if (axisList[i].class != "R") && (proposals[j].class != "R")
                    proposals[j].degree = lcm(axisList[i].degree, proposals[j].degree);
                    proposals[j].degree = min(proposals[j].degree, 20)
                end
                if (axisList[i].class != proposals[j].class)
                    if proposals[j].class=="R"
                        proposals[j].coordinate = axisList[i].coordinate
                    end
                    proposals[j].class = "D"
                end
            elseif dist[1] < 0.3
                if (axisList[i].degree > proposals[j].degree) && (proposals[j].degree <= 2)
                    temp = deepcopy(proposals[j])
                    proposals[j] = deepcopy(axisList[i])
                    axisList[i] = temp
                end
                if axisList[i].degree==2
                    proposals[j].class="D"
                elseif axisList[i].degree==-1 
                    if proposals[j].degree==-1
                        proposals[j].coordinate[:,2]=cross(vec1,vec2)
                        proposals[j].coordinate[:,2]=proposals[j].coordinate[:,2] ./ norm(proposals[j].coordinate[:,2])
                        proposals[j].coordinate[:,3]=vec1 ./ norm(vec1)
                        proposals[j].coordinate[:,1]=cross(proposals[j].coordinate[:,2], proposals[j].coordinate[:,3])
                        proposals[j]=Axis("C", 2, proposals[j].coordinate)
                    else
                        proposals[j].coordinate[:,3]=cross(vec1,vec2)
                        proposals[j].coordinate[:,3]=proposals[j].coordinate[:,3] ./ norm(proposals[j].coordinate[:,3])
                        proposals[j].coordinate[:,1]=cross(proposals[j].coordinate[:,2], proposals[j].coordinate[:,3])
                    end
                else
                    merged = false
                end
            end
        end
        if !merged
            push!(proposals, axisList[i])
        end
    end
    return proposals
end

for idx = 1:1  #size(models,1)-1
#    md5 = models[idx];
    md5 = ARGS[2]
    println(md5);
    #try
        p = sortperm(vec(NNgraph[idx,:]));
        sym_anchor = readSymmetry(synsetID, md5);
        axisList = [];
        if (p[end] > 0)
            for n = 2:16
                md5_1 = models[p[n]];
                #try
                    sym_neighbor = readSymmetry(synsetID, md5_1);
                    #println(md5_1)
                    #println(sym_neighbor)
                    append!(axisList, sym_neighbor);
                #catch err1
                #    if !isa(err1, LoadError)
                #        throw(err1)
                #    end
                #end
            end
        end
        if !isempty(axisList)
            PC = loadPC(synsetID, md5);
            proposals= mergeAxis(axisList, sym_anchor);
            println(axisList)
            println(sym_anchor)
            println(proposals)
            logfile = open("Results/"*synsetID*"/"*md5*".log2", "w");
            if size(proposals,1)==1
                symtype, translate, coordinate = refineSym(proposals[1], proposals[1].coordinate, PC, logfile);
            else
                num_high_deg = 0
                highest_deg = 1
                highest_coordinate = proposals[1].coordinate
                highest_symtype = "E"
                highest_translate = zeros(3,1)
                for i = 1:size(proposals,1)
                    symtype, translate, coordinate = refineSym(proposals[i], proposals[i].coordinate, PC, logfile)
                    println(symtype)
                    println(coordinate)
                    degree = tryparse(Int64, symtype[2:end]);
                    if isnull(degree)
                        degree = tryparse(Int64, symtype[2:end-1]);
                    end
                    if isnull(degree)
                        degree = 1
                    else
                        degree = get(degree);
                    end
                    if degree > 2
                        num_high_deg += 1
                    end
                    if degree > highest_deg
                        highest_deg = degree
                        highest_coordinate = coordinate
                        highest_symtype = symtype
                    end
                    highest_translate += translate
                    PC = PC .+ translate'
                end
                if num_high_deg > 1
                    if highest_deg == 4
                        symtype = "O"
                    elseif highest_deg == 5
                        symtype = "I"
                    elseif highest_deg == 3
                        symtype = "T"
                    else
                        symtype = "O3"
                    end
                else
                    symtype = highest_symtype
                end
                translate = highest_translate
                coordinate = highest_coordinate
            end
            close(logfile)
            file = "Results/"*synsetID*"/"*md5*".sym3t_g"
            saveSymmetry(file, symtype, translate, coordinate)
        end

    #catch err
    #    log = open(synsetID*".errorlog", "a");
    #    @printf(log, "%s\n", md5);
    #    @printf(log, "%s\n", err);
    #    close(log);
    #end

end





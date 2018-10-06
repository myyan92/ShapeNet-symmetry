using IOUtil
using RefineAxisLib

function generateAxis(symType, translate, coordinate)
    if symType=="E" 
        return []; 
    end
    if symType=="Cs"
        return [Axis("R", -1, coordinate)];
    end
    if (symType[1]=='C' || symType[1]=='D')
        degree = tryparse(Int64, symType[2:end]);
        if isnull(degree)
            degree = parse(Int64, symType[2:end-1]);
        else 
            degree = get(degree);
        end
    end
    if symType[1]=='C'
        if symType[end]=='h'
            return [Axis("D", degree, coordinate)];
        else
            return [Axis("C", degree, coordinate)];
        end
    end
    if symType[1]=='D'
        return [Axis("D", degree, coordinate)];
    end
    if symType[1]=='O'
        if symType=="O3"
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
    error("symmetry not implemented")    
end


function mergeAxis(nbAxises, selfAxises)
    if isempty(nbAxises)
        return selfAxises
    end
    if !isempty(selfAxises)
        proposals = deepcopy(selfAxises);
    else
        proposals = [deepcopy(nbAxises[1])]
    end
    for i = 1:size(nbAxises,1)
        merged = false
        if nbAxises[i].class=="E"
            merged = true
        end
        for j = 1:size(proposals,1)
            if merged continue; end
            merged = true
            if proposals[j].class=="E"
                proposals[j]=deepcopy(nbAxises[i])
                continue
            end
            vec1 = proposals[j].coordinate[:,2]
            if proposals[j].class=="R"
                vec1 = proposals[j].coordinate[:,1]
            end
            vec2 = nbAxises[i].coordinate[:,2]
            if nbAxises[i].class=="R"
                vec2 = nbAxises[i].coordinate[:,1]
            end
            dist = abs(vec2'*vec1);
            if dist[1] > 0.97
                # RR: duplicate, CC: merge degree, CR / RC: becomes Cnh group (D verify), copy coordinate from the C axis
                if (nbAxises[i].class != "R") && (proposals[j].class != "R")
                    proposals[j].degree = lcm(nbAxises[i].degree, proposals[j].degree);
                end
                if (nbAxises[i].class != proposals[j].class)
                    if proposals[j].class=="R"
                        proposals[j] = deepcopy(nbAxises[i])
                    end
                    proposals[j].class = "D"
                end
            elseif (nbAxises[i].class == "R") && (proposals[j].class != "R") 
                # from C to Cv, or D to Dh
                prod = abs(vec2'*proposals[j].coordinate[:,1]);
                clamp!(prod, 0, 1)
                angle_err = acos(prod[1]) / 3.14159 * 180
                while (angle_err > 90.0 / proposals[j].degree)  angle_err = angle_err - 180.0 / proposals[j].degree end
                if (abs(angle_err)>10) merged=false end
            elseif (nbAxises[i].class != "R") && (proposals[j].class == "R") 
                # from C to Cv, or D to Dh
                prod = abs(vec1'*nbAxises[i].coordinate[:,1]);
                clamp!(prod, 0, 1)
                angle_err = acos(prod[1]) / 3.14159 * 180
                while (angle_err > 90.0 / nbAxises[i].degree)  angle_err = angle_err - 180.0 / nbAxises[i].degree  end
                if (abs(angle_err)>10) 
                    merged=false
                else
                    proposals[j] = deepcopy(nbAxises[i])
                end
            elseif (nbAxises[i].degree == 2) && (proposals[j].class != "R") 
                # from C/D to D
                prod = abs(vec2'*proposals[j].coordinate[:,3]);
                clamp!(prod, 0, 1)
                angle_err = acos(prod[1]) / 3.14159 * 180
                while (angle_err > 90.0 / proposals[j].degree)  angle_err = angle_err - 180.0 / proposals[j].degree end
                if (abs(angle_err)>10)
                    merged=false
                else
                    proposals[j].class = "D"
                end
            elseif (nbAxises[i].class != "R") && (proposals[j].degree == 2)
                # from C/D to D
                prod = abs(vec1'*nbAxises[i].coordinate[:,3]);
                clamp!(prod, 0, 1)
                angle_err = acos(prod[1]) / 3.14159 * 180
                while (angle_err > 90.0 / nbAxises[i].degree)  angle_err = angle_err - 180.0 / nbAxises[i].degree end
                if (abs(angle_err)>10)
                    merged=false
                else
                    proposals[j] = deepcopy(nbAxises[i])
                end
            end
            if !merged
                # println(proposals[j])
                # println(nbAxises[i])
            end  
        end
        if !merged
            push!(proposals, nbAxises[i])
        end
    end
    return proposals
end



synsetID = ARGS[1];
model_list = "./deduplicate_lists/"*synsetID*".txt";
models = readall(model_list);
models = split(models, '\n');
size(models,1)

graphFile = "./NNGraph/"*synsetID*".txt";
NNgraph = readdlm(graphFile);
size(NNgraph,1)


for idx = 1:size(models,1)-1
    md5 = models[idx];
    # println(md5);
    p = sortperm(vec(NNgraph[idx,:]));
    symType, translate, coordinate = readSymmetry(synsetID, md5);
    self_sym = generateAxis(symType, translate, coordinate)
    neighborList = [];
    for n = 2:16
        md5_1 = models[p[n]];
        symType_n, translate_n, coordinate_n = readSymmetry(synsetID, md5_1);
        sym_neighbor = generateAxis(symType_n, translate_n, coordinate_n)
        append!(neighborList, sym_neighbor);
    end
    proposals= mergeAxis(neighborList, self_sym);
    if !isempty(proposals)
        Mesh = loadMesh(synsetID, md5);
        Mesh.vertices = Mesh.vertices .- translate'
        densityProposal = GetSampleDensityProposal(Mesh, 8000);
        PC = SamplePoints(Mesh, densityProposal);
        logfile = open("Results/"*synsetID*"/"*md5*".log_g", "w");

        for i = 1:3  # refine and merge multiple times to be accurate
            for i = 1:size(proposals,1)
                symtype, degree, translate, coordinate = refineSym(proposals[i], proposals[i].coordinate, PC, logfile);
                proposals[i].class = symtype
                if proposals[i].class == "Cs"
                    proposals[i].class = "R"
                elseif proposals[i].class[1]=='C' && proposals[i].class[end]=='h'
                    proposals[i].class = "D"
                else
                    proposals[i].class = proposals[i].class[1:1]
                end
                proposals[i].degree = degree
                proposals[i].coordinate = coordinate
            end
            proposals= mergeAxis(proposals, [])
        end
        close(logfile)
    end
    if size(proposals, 1)>1
        println(md5)
    end
    for i = 1:size(proposals,1)
        symtype, degree, translate, coordinate = refineSym(proposals[i], proposals[i].coordinate, PC, logfile);
        file = "Results/"*synsetID*"/"*md5*".sym3t_g$(i)"
        saveSymmetry(file, symtype, translate, coordinate)
    end         

end





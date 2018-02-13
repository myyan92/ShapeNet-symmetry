using MeshIO
using FileIO
using Distances
using NearestNeighbors
@everywhere include("icp.jl")
@everywhere include("symmetrySpectral.jl")
dataroot = "../../symmetrySpectral/results/Lamp/"
filelist = readdir(dataroot);
isoff(x) = contains(x, ".off")
filter!(isoff, filelist)
filelist = [dataroot * f[1:end-4] for f in filelist]
pmap(symmetrySpectral, filelist)

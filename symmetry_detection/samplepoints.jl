using FileIO
using MeshIO

function ComputeFaceArea(p1, p2, p3)
    @assert length(p1) == 3 && length(p2) == 3 && length(p3) == 3    
    p1 = vec(p1)
    p2 = vec(p2)
    p3 = vec(p3)
    A = cross(p2 - p1, p3 - p1)
    area = 0.5 * norm(A)
end

function ComputeMeshArea(vertices, faces)
    area = 0
    for f in faces
        pt1 = [vertices[f[1]+1][1], vertices[f[1]+1][2], vertices[f[1]+1][3]]
        pt2 = [vertices[f[2]+1][1], vertices[f[2]+1][2], vertices[f[2]+1][3]]
        pt3 = [vertices[f[3]+1][1], vertices[f[3]+1][2], vertices[f[3]+1][3]]
        A = cross(pt2 - pt1, pt3 - pt1)
        area += 0.5 * norm(A)
    end
    area
end

function GetSampleDensityProposal(vertices, faces, expectedNumPoints)
    expectedNumPoints / ComputeMeshArea(vertices, faces)
end

function SamplePoints(vertices, faces, numPoints)
    densityProposal = GetSampleDensityProposal(vertices, faces, numPoints)
    if isinf(densityProposal)        
        return vertices
    end
    
    points = []
    for f in faces
        p1 = [vertices[f[1]+1][1], vertices[f[1]+1][2], vertices[f[1]+1][3]]
        p2 = [vertices[f[2]+1][1], vertices[f[2]+1][2], vertices[f[2]+1][3]]
        p3 = [vertices[f[3]+1][1], vertices[f[3]+1][2], vertices[f[3]+1][3]]
        
        area = ComputeFaceArea(p1, p2, p3)
        expectedSampleNum = area * densityProposal        
        actualSampleNum = floor(Int, expectedSampleNum)
        epsilon = rand()
        if epsilon <= expectedSampleNum - actualSampleNum
            actualSampleNum += 1
        end
        
        # Uniform sampling in a triangle. The key idea is to generate a uniformly
        # distributed triplet (alpha, beta, gamma) in a triangle defined by:
        # 0 <= alpha <= 1, 0 <= beta <= 1, 0 <= gamma <= 1.
        # alpha + beta + gamma = 1.
        # Equivalently, we generate (alpha, beta) that are uniformly distributed
        # in a triangle defined by:
        # 0 <= alpha <= 1, 0 <= beta <= 1, alpha + beta <= 1.
        for j = 1:actualSampleNum            
            # Generate uniform (alpha, beta).
            alpha = rand()
            beta = rand()
            if (alpha + beta > 1) 
                alpha = 1 - alpha
                beta = 1 - beta
            end
            # Generate uniform (alpha, beta, gamma).
            gamma = 1 - alpha - beta

            # Barycentric coordinates.
            newPoint = p1 * alpha + p2 * beta + p3 * gamma
            push!(points, newPoint)            
        end
    end
    retPoints = zeros(length(points), 3)
    for i = 1:length(points)
        retPoints[i, :] = points[i]
    end
    retPoints
end

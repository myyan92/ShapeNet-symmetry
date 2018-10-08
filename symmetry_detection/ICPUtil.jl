module ICPUtil

using Distances
using NearestNeighbors

export icp

function match_bruteForce(q, p)
    m = size(p,2);
    n = size(q,2);    
    match = zeros(1,m);
    mindist = zeros(1,m);
    dists = pairwise(Euclidean(), q, p)
    mindist, match = findmin(dists, 1)
    match = (match-1) .% n + 1
    match, mindist
end

function match_kDtree(q, p, kdtree)
    match, mindist = knn(kdtree,p,1);
    match = [x[1] for x in match]
    mindist = [x[1] for x in mindist]
    match, mindist
end

function eq_point(q,p)
    m = size(p,2);
    n = size(q,2);

    # find data centroid and deviations from centroid
    q_bar = sum(q, 2) / n;
    q_mark = broadcast(-, q, q_bar);
    # find data centroid and deviations from centroid
    p_bar = sum(p, 2) / m;
    p_mark = broadcast(-, p, p_bar);
    
    N = p_mark*q_mark'; 
    U,S,V = svd(N); 
    R = V*diagm([1, 1, det(U*V')])*U';
    T = q_bar - R*p_bar;
    
    R, T
end

function eq_plane(q,p,n)

    c = zeros(3, size(q,2))
    for i = 1:size(q,2)
        c[:,i] = cross(vec(p[:,i]),vec(n[:,i]));
    end

    cn = vcat(c,n);
    C = cn*cn';

    b = -[sum((p-q).*repmat(cn[1,:],3,1).*n);
          sum((p-q).*repmat(cn[2,:],3,1).*n);
          sum((p-q).*repmat(cn[3,:],3,1).*n);
          sum((p-q).*repmat(cn[4,:],3,1).*n);
          sum((p-q).*repmat(cn[5,:],3,1).*n);
          sum((p-q).*repmat(cn[6,:],3,1).*n)];

    X = C\b;
    cx = cos(X[1]); cy = cos(X[2]); cz = cos(X[3]); 
    sx = sin(X[1]); sy = sin(X[2]); sz = sin(X[3]); 

    R = [cy*cz cz*sx*sy-cx*sz cx*cz*sy+sx*sz;
         cy*sz cx*cz+sx*sy*sz cx*sy*sz-cz*sx;
         -sy cy*sx cx*cy];
    T = X[4:6];
    
    R, T
end

function rms_error(p1,p2)
    dsq = sum((p1 - p2).^2,1);
    ER = sqrt(mean(dsq));
    ER
end

function lsqnormest(p, k)
    m = size(p,2)
    n = zeros(3,m)

    kdtree = KDTree(p)
    neighbors, dists = knn(kdtree, p, k+1, true)

    for i = 1:m
        x = p[:,neighbors[i][2:end]]
        p_bar = 1/k * sum(x, 2)
        x = broadcast(-, x, p_bar)
        P = x * x';

        D,V = eig(P)
        idx = indmin(D) 
        n[:,i] = V[:,idx];   
    end
    n
end

function icp(q,p,k=10, normal=[]; Matching="kDtree", Minimize="plane", ReturnAll=false, WorstReject=0)
"""
     Perform the Iterative Closest Point algorithm on three dimensional point
  clouds.

  [TR, TT, ER] = icp(q,p,k)   returns the rotation matrix TR and translation
  vector TT that minimizes the distances from (TR * p + TT) to q. p is a 3xm matrix
     and q is a 3xn matrix. k is the number of iterations. ER is the RMS of errors for k
  iterations in a (k+1)x1 vector. ER(0) is the initial error. Also returns the calculation times per
  iteration in a (k+1)x1 vector.

  Additional settings:

  Matching
      {bruteForce} | kDtree
      Specifies how point matching should be done. 

  Minimize
      point | {plane}
      Defines whether point to point or point to plane minimization
      should be performed. 

  Normals
      A matrix of normals for the n points in q.

  ReturnAll
      {false} | true
      Determines whether R and T should be returned for all iterations
      or only for the last one. If this option is set to true, R will be
      a 3x3x(k+1) matrix and T will be a 3x1x(k+1) matrix.

  WorstReject
      {0} | scalar in [0; 1]
      Reject worst point pairs, based on their Euclidean distance > threshold.

  Martin Kjer and Jakob Wilm, Technical University of Denmark, 2012

"""
    # Allocate vector for time and RMS of errors in every iteration.
    ER = zeros(k+1,1);

    Np = size(p,2);

    # Transformed data point cloud
    pt = p;

    # Initialize temporary transform vector and matrix.
    T = zeros(3,1);
    R = eye(3,3);

    # Initialize total transform vector(s) and rotation matric(es).
    TT = zeros(3,1, k+1);
    TR = zeros(3,3, k+1);
    TR[:,:,1] = eye(3,3);

    # If Minimize == "plane", normals are needed
    if (Minimize == "plane") && isempty(normal)
        normal = lsqnormest(q,20);
    end

    # If Matching == 'kDtree', a kD tree should be built (req. Stat. TB >= 7.3)
    if Matching == "kDtree"
        kdtree = KDTree(q);
    end

    # Go into main iteration loop
    for dk=1:k
        if Matching == "kDtree"
            match, mindist = match_kDtree(q,pt,kdtree); 
        else 
            match, mindist = match_bruteForce(q,pt);
        end
        p_idx = trues(1, Np);
        q_idx = match;
        # If worst matches should be rejected
        if WorstReject>0
            edge = round(Int, (1-WorstReject)*sum(p_idx));
            idx = sortperm(mindist);
            for i in idx[edge:end]
                p_idx[i] = false;
            end
            q_idx = match[vec(p_idx)];
            mindist = mindist[vec(p_idx)];
        end

        if dk == 1
            ER[dk] = sqrt(sum(mindist.^2)/length(mindist));
        end

        if Minimize == "point"
            R,T = eq_point(q[:,vec(q_idx)],pt[:,vec(p_idx)]);
        elseif Minimize == "plane"
            R,T = eq_plane(q[:,vec(q_idx)],pt[:,vec(p_idx)],normal[:,vec(q_idx)]);
        end

        # Add to the total transformation
        TR[:,:,dk+1] = R*TR[:,:,dk];
        TT[:,:,dk+1] = R*TT[:,:,dk]+T;

        # Apply last transformation
        pt = TR[:,:,dk+1] * p + repmat(TT[:,:,dk+1], 1, Np);
        ER[dk+1] = rms_error(q[:,vec(q_idx)], pt[:,vec(p_idx)]);
        if ER[dk]-ER[dk+1] < 1e-4 
            TR = TR[:,:,1:dk+1]
            TT = TT[:,:,1:dk+1]
            ER = ER[1:dk+1]
            break; 
        end
    end

    if !(ReturnAll)
        TR = TR[:,:,end];
        TT = TT[:,:,end];
    end

    TR, vec(TT), ER
    
end



end # module end

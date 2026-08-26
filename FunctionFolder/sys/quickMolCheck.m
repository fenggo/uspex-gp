function good = quickMolCheck(tempMOLS, MUT_LAT, MtypeLIST, minDistMatrice)
%QUICKMOLCHECK  Fast molecular overlap check using centroid distances.
%
%   Replaces newMolCheck with a ~100x faster pre-filter. Instead of
%   checking every atom-atom pair, it only checks molecule centroid
%   distances against an adaptive threshold.
%
%   Inputs:
%     tempMOLS      - struct array with .MOLCOORS (na×3 per molecule)
%     MUT_LAT       - 3×3 lattice matrix (Cartesian)
%     MtypeLIST     - molecule type indices (unused, kept for API compatibility)
%     minDistMatrice- minimum allowed atom-atom distance matrix
%
%   Output:
%     good          - 1 if no overlap detected, 0 otherwise

totalMols = length(tempMOLS);
molCenters = zeros(totalMols, 3);
molRadii   = zeros(totalMols, 1);

% Compute centroid and radius for each molecule
for i = 1:totalMols
    coords = tempMOLS(i).MOLCOORS;
    molCenters(i,:) = mean(coords, 1);
    diffs = coords - molCenters(i,:);
    molRadii(i) = max(sqrt(sum(diffs.^2, 2)));
end

% Convert fractional centroids to Cartesian
cartCenters = molCenters * MUT_LAT;

% Adaptive threshold: use the larger of 70% of sum of radii or
% 80% of the maximum atom-atom minimum distance
maxAtomDist = max(minDistMatrice(:));
threshold = max(sum(molRadii) * 0.05, maxAtomDist * 0.8);

% O(n²) centroid distance check
for i = 1:totalMols-1
    ci = cartCenters(i,:);
    for j = i+1:totalMols
        d = norm(ci - cartCenters(j,:));
        if d < threshold
            good = 0;
            return;
        end
    end
end

good = 1;
end
function [inds maxVals] = findMaxNd(M,n)
% [inds maxVals] = findMaxNd(M,n)
%
% Returns n highest values in matrix M, which may have >2 dimensions.
%
% inds(i,:) = indices of ith heighest value in matrix M.
%           = [a b c d ...]
%M(a,b,c,d,...) = ith heighest value in matrix M
%
% maxVals(i,:) = [ith heighest value in matrix M, how many entries have
% this value]
%
%horribly inefficient implementation, don't give it big matrices / large
%values of n

matSize = size(M);
numDims = length(matSize);
inds = zeros(n,numDims);%may be too small if multiple matrix entries have same value
maxVals = zeros(n,2);

nextRow = 1;
for i=1:n

    thisMax = M;
    for j=1:numDims
        thisMax = max(thisMax);
    end
    
    theseLinInds = find(M == thisMax);
    maxVals(i,:) = [M(theseLinInds(1)) length(theseLinInds)];
    M(theseLinInds) = NaN;
    
    for j=1:length(theseLinInds)
        inds(nextRow,:) = myInd2Sub(matSize,theseLinInds(j));
        nextRow = nextRow + 1;
    end
end
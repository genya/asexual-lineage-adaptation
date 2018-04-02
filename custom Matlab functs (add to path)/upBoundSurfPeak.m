function upperPeak = upBoundSurfPeak(surf,paramRange)
%upperPeak = upBoundSurfPeak(surf,paramRange)
%
% surf = N-dimensional matrix
%
% paramRange = cell array length N
% paramRange{i} = parameter range of dimension i for matrix surf
%
% The peak of surf needs to be at least 2 points away from the edge along 
% all axes or else "Index exceeds matrix dimensions" errors will occur.
%
% surf may only have a single entry equal to the max value for the matrix.
%
% surf may not contain NaNs in the vicinity of peak.

%find the peak 
[peakInds, ~] = findMaxNd(surf,1);

if size(peakInds,1) > 1
    error('multiple peaks')
end
peakLinInd = mySub2Ind(size(surf),peakInds);

%compute derivatives
numDims = length(paramRange);
surfSize = size(surf);

firDivs = zeros(numDims,1);
secDivs = zeros(numDims,numDims);
% firDivsPlus = zeros(numDims,numDims);
% firDivsMinus = zeros(numDims,numDims);
    %firDivs[PlusMinus](i,j) = 
    %   first derivative wrt param i at the pt neighboring the peak in the
    %   [plus/minus] direction along the jth axis.

for i=1:numDims
    
    firDivs(i) = firDivFunct(surf,peakInds,i);   
    
    for j=1:numDims
        neighbIndPlus = peakInds;
        neighbIndPlus(j) = neighbIndPlus(j) + 1;
        firDivsPlus = firDivFunct(surf,neighbIndPlus,i);
        
        neighbIndMinus = peakInds;
        neighbIndMinus(j) = neighbIndMinus(j) - 1;
        firDivsMinus = firDivFunct(surf,neighbIndMinus,i);        
               
        secDivs(i,j) = (firDivsPlus - firDivsMinus)...
            /( paramRange{j}(neighbIndPlus(j)) - paramRange{j}(neighbIndMinus(j)));
        
    end
end


%define paraboloid
thisPoly = @(x) -( ...
    surf(peakLinInd) + sum(x(:).*firDivs)...
    + (1/2)*sum(sum( xMat(x(:)).*secDivs )) );

    %should I check that this paraboloid is concave down? 

%find max
options = optimset('fminsearch');
%  options = optimset(options,'MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',10^-3,'TolX',10^-4);
[~, upperPeak] = fminsearch(thisPoly,zeros(numDims,1),options);
upperPeak = -upperPeak;

%sub-functions
    function prodMat = xMat(x)
        prodMat = zeros(length(x),length(x));
        thisX = x(:);
        for xMatInd=1:length(x)
            prodMat(:,xMatInd) = x(xMatInd)*thisX;
        end
    end

    function firstDeriv = firDivFunct(thisSurf,ptInds,dimInd)
        peakIndPlus = ptInds;
        peakIndPlus(dimInd) = ptInds(dimInd) + 1;
        peakIndMinus = ptInds;
        peakIndMinus(dimInd) = ptInds(dimInd) - 1;
        
        fPlus = thisSurf(mySub2Ind(surfSize,peakIndPlus));
        fMinus = thisSurf(mySub2Ind(surfSize,peakIndMinus));

        firstDeriv = (fPlus - fMinus)...
            /( paramRange{dimInd}(peakIndPlus(dimInd)) - paramRange{dimInd}(peakIndMinus(dimInd)) );
    end
        
end
function linearInd = mySub2Ind(matSize, subscripts)
%linearInd = mySub2Ind(matSize, subscripts)
%
% returns the same value as IND = sub2ind(SIZ,I1,I2,...,IN)
% but for input [I1,I2,...,IN] where each entry is a scalar, i.e.
% 
% linearInd = mySub2Ind(SIZ, [I1,I2,...,IN] ) = sub2ind(SIZ,I1,I2,...,IN)
% 

subscripts = subscripts(:);
matSize = matSize(:);

if length(subscripts) == 1
    linearInd = subscripts;
else
    dimWeights = cumprod(matSize);
    dimWeights = dimWeights(1:end-1);
    linearInd = subscripts(1) + sum((subscripts(2:end)-1).*dimWeights(:));
end
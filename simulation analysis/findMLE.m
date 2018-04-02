function [ML, ML_ind] = findMLE(trajTypeCnts,trajPeakSdown,trajSimTot,trajInd)
%function [ML, ML_ind] = findMLE(trajTypeCnts,trajPeakSdown,trajSimTot)

if nargin==3
    trajInd = [1:length(trajPeakSdown)]';
end

logLikeList = zeros(size(trajTypeCnts));

%likelihood estimate & uncertainty based on Wilson Score
zScore2 = 2^2; % roughly alpha=5%, see 2*(1-normcdf(zScore,0,1))
likEst = @(sawDat,simsSoFar) (sawDat + zScore2/2)./(simsSoFar + zScore2);

%iterate through trajTypeCnts

pick = sub2ind( size(trajTypeCnts{1}),...
        trajPeakSdown(:),...
        trajInd(:));
for i=1:length(trajTypeCnts)
    %compute log likelihood sum 
    logLikeList(i) = sum(log(...
        likEst(trajTypeCnts{i}(pick), trajSimTot{i}(:))...
        ));
end

% where's the MLE?
[ML, ML_ind] = max(logLikeList);
function [pVal logLikeDist passQC] = pValFromMLE(trials,dataMLE_height,dataMLE_ind,trajTypeDist,trajTypeCnts,trajSimTot,datSize,datSize2,peakRegInteriorInds,peakRegBorderInds)
%[pVal logLikeDist passQC] = pValFromMLE(trials,dataMLE_height,dataMLE_ind,trajTypeDist,trajTypeCnts,trajSimTot,datSize,datSize2,peakRegInteriorInds,peakRegBorderInds)
%
% computes p-value using the likelihood as a test statistic. 
%  = fraction 
%
% logLikeDist = a list of maximum likelihoods of length trials, each of 
% which is the maximum likelihood of synthetic data generated from
% trajTypeDist{dataMLE_ind} over the space of all of trajTypeDist.
%
% passQC(i) = 1 => 
%  a) MLE of logLikeDist(i) was in the interior of the peak region (least 2 cells from border) 
%  b) its uncertainty was less than 1/2
%  c) border of peak region is at least 5 log likelihood units below ML
%
% notes on inputs:
% peakRegInteriorInds = indices of points in peak reg interior
% peakRegBorderInds = indices of peak region border
% these are logical indices
 
%seed random number generator
rng('shuffle');

logLikeDist = zeros(trials,1);
numChecks = 3;
passQC = false(numChecks,trials);
    
%likelihood estimate & uncertainty based on Wilson Score
zScore2 = 4; % roughly alpha=5%, see 2*(1-normcdf(zScore,0,1))
likEst = @(sawDat,simsSoFar) (sawDat + zScore2/2)./(simsSoFar + zScore2);
likSig = @(sawDat,simsSoFar) sqrt( sawDat.*(1-sawDat./simsSoFar) + zScore2/4 )./(simsSoFar + zScore2);


synthPeakSdownBins = zeros(size(trajTypeDist{1},2),1);
synthLogLikes = zeros(length(trajTypeDist),1);
synthLogLikeSigs = zeros(length(trajTypeDist),1);

for j=1:trials
    
    %generate synthetic data at MLE given the experimental data
    for i=1:length(synthPeakSdownBins)
        synthPeakSdownBins(i) = randBeans(trajTypeDist{dataMLE_ind}(:,i),1);
    end
    
    %for the synthetic data set, find the MLE
    pick = sub2ind(datSize,synthPeakSdownBins,datSize2);
    for i=1:length(trajTypeDist)        
        synthLogLikes(i) =  sum(log(...
           likEst(trajTypeCnts{i}(pick),trajSimTot{i}(:))...
           ));  
       synthLogLikeSigs(i) = sqrt(sum((...
           likSig(trajTypeCnts{i}(pick),trajSimTot{i}(:))...
           ./synthLogLikes(i)...
           ).^2));
    end  
    
    [logLikeDist(j) maxSynthLogLikInd] = max(synthLogLikes);  
    
    
    thisQC = false(numChecks,1);
    %check that its uncertainty is less than <1/2 log likelihood units
     thisQC(1) = synthLogLikeSigs(maxSynthLogLikInd) < 1/2;
    
    %check that it is in the peak region
    thisQC(2) = peakRegInteriorInds(maxSynthLogLikInd) == 1;
    
    %check that boundaries of the peak region have log likelihood at least
    %2 units lower
    thisQC(3) = sum(logLikeDist(j) - 2 < synthLogLikes(peakRegBorderInds))==0;
    
    passQC(:,j) = thisQC(:);
          
end

%compute p-value = fraction of max likelihoods less than ML_data
pVal = sum( logLikeDist < dataMLE_height)/trials;
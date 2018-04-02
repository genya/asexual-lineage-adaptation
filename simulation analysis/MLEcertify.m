%% load sim results [mandatory to run this cell]
clear
homeDir = 'C:\Users\genya\Google Drive\BMA sim code\for distribution\';
cd(homeDir);

%specify what simulations to load 
% simBatchInd = 23;
% DFE_form = 'unif';

% simBatchInd = 14;
% DFE_form = 'exp';

simBatchInd = 17;
DFE_form = 'dir';

% simBatchInd = 6;
% DFE_form = 'truncExp';

loadThese = [DFE_form 'DFEresults_' num2str(simBatchInd)];
load(['simulation results\' ...
	loadThese])

%% check for results with param differences below significance [mandatory to run this cell]
%(spurious differences may be introduced by Matlab bull shit, e.g. "guard digits")

combos = zeros(size(param2results,1),1);
numSigFigs = 5;

comboNumber =1;
for i=1:size(param2results,1)
 
   pick = true(size(param2results,1),1); 
   for j=1:size(param2results,2)
        pick = pick & ...
            abs(param2results(:,j) - param2results(i,j)) <= 10^-numSigFigs;
   end
   
   if sum(pick) > 1
       combos(pick)  = comboNumber;
       comboNumber = comboNumber + 1;
   end
end
if isempty(find(combos, 1))
    clear pick* i combo* j
else
    display('need to combine redundant sim results')
    keyboard
end

for i=1:length(DFEparamRange)
    for j=1:length(DFEparamRange{i})-1
        for k=j+1:length(DFEparamRange{i})
            if abs(DFEparamRange{i}(j)-DFEparamRange{i}(k)) < 10^-numSigFigs
                DFEparamRange{i}(k) = DFEparamRange{i}(j);
            end
        end
    end
    DFEparamRange{i} = unique(DFEparamRange{i});
end
clear i j k pick
%% compute log likelihoods [mandatory to run this cell]

trajFeature = 1; %if =1, use peakFreqSdown, if =0 use sweep or not;

if trajFeature == 0
    trajSweepOrNot = zeros(size(trajPeakSdown));
    trajSweepOrNot(trajPeakSdown == 2) = 1; 
    trajSweepOrNot(trajSweepOrNot == 0) = 2;
end

logLikeList = zeros(size(trajTypeCnts));
logLikeSigList = zeros(size(trajTypeCnts));
trajTypeDist = cell(size(trajTypeCnts));
trajSimTot = cell(size(trajTypeCnts));

%likelihood estimate & uncertainty based on Wilson Score
zScore2 = 2^2; % roughly alpha=5%, see 2*(1-normcdf(zScore,0,1))
likEst = @(sawDat,simsSoFar) (sawDat + zScore2/2)./(simsSoFar + zScore2);
likSig = @(sawDat,simsSoFar) sqrt( sawDat.*(1-sawDat./simsSoFar) + zScore2/4 )./(simsSoFar + zScore2);

%iterate through trajTypeCnts
oneTwoThree = [1:length(trajPeakSdown)]';

for i=1:length(trajTypeCnts)

    %do all dat trajs have some simulations?
    theseSimTotals = sum(trajTypeCnts{i});    
    trajSimTot{i} = theseSimTotals(:);
    theseSimTotals = theseSimTotals(:);
   
    if find(theseSimTotals ==0,1)
        logLikeList(i) = NaN;
        continue;
    end

    %compute likelihood and uncertainty
    if trajFeature==1
        pick = sub2ind( size(trajTypeCnts{i}),...
            trajPeakSdown(:),...
            oneTwoThree(:));
        theseLikes = likEst(trajTypeCnts{i}(pick), theseSimTotals(oneTwoThree));
        theseSigs = likSig(trajTypeCnts{i}(pick), theseSimTotals(oneTwoThree));     
    elseif trajFeature==0
        theseTrajTypeCnts = [trajTypeCnts{i}(2,oneTwoThree); theseSimTotals(oneTwoThree)'-trajTypeCnts{i}(2,oneTwoThree)];%assumes peakFreqSdownBins(2,:) = [1 NaN];
        pick = sub2ind( size(theseTrajTypeCnts),...
            trajSweepOrNot(:),...
            oneTwoThree(:));
        theseLikes = likEst(theseTrajTypeCnts(pick), theseSimTotals(oneTwoThree));
        theseSigs = likSig(theseTrajTypeCnts(pick), theseSimTotals(oneTwoThree));            
    else error('trajFeature undefined')
    end
    
    %sum log likihood & propagate uncertainty
    logLikeList(i) = sum(log(theseLikes));
    logLikeSigList(i) = sqrt(sum( (theseSigs./theseLikes).^2));
       
    %normalize traj type cnts to get distribution of traj types
    trajTypeDist{i} = zeros(size(trajTypeCnts{i}));
    for j=1:size(trajTypeCnts{i},2)
        trajTypeDist{i}(:,j) = trajTypeCnts{i}(:,j)/sum(trajTypeCnts{i}(:,j));
    end
   
end

% where's the MLE?
[ML, MLE_ind] = max(logLikeList);
ML
MLE_params = param2results(MLE_ind,:)
peakInds = zeros(size(DFEparamRange));
for i=1:length(peakInds)
    peakInds(i) = find(...
        abs(DFEparamRange{i} - param2results(MLE_ind,i)) <= 10^-numSigFigs...
        );
end

clear numDims oneTwoThree sizTrajTypeCnts
clear this* these* pick i j zScore* likEst likSig k
%% draw likelihood surface [mandatory to run this cell]
 
% initialize surface variables
sizeSurf = zeros(1,length(DFEparamRange));
for i=1:length(DFEparamRange)
    sizeSurf(i) = length(DFEparamRange{i});
end
clear i

likeSurf = nan(sizeSurf);
likeSigs = nan(sizeSurf);
simTots = zeros(sizeSurf);

%iterate through logLikeList
for i=1:length(logLikeList)
       
    %find location on likelihood surface
    theseInds = zeros(size(param2results(i,:)));
    for j=1:length(theseInds)
       theseInds(j) =  find(...
           abs( DFEparamRange{j} - param2results(i,j) ) <= 10^-numSigFigs ...
            );
    end
    thisInd = mySub2Ind(size(likeSurf),theseInds);

    if find( trajSimTot{i} < 1,1)
        likeSurf(thisInd) = NaN;
        likeSigs(thisInd) = NaN;  
        continue;
    end    
    
    likeSurf(thisInd) = logLikeList(i);
    likeSigs(thisInd) = logLikeSigList(i);
    simTots(thisInd) = sum(trajSimTot{i});    
end
clear i j this* these*



%% eyeball likelihood surface [optional to run this cell]

plotThis1 = likeSurf(:,:);
plotThis2 = likeSigs(:,:);
plotThis3 = log10(simTots(:,:)...
    *length(sInoc)/length(trajs));
xVar = 2; yVar = 1; 
plotX = DFEparamRange{xVar};
plotY = DFEparamRange{yVar};

figure;
subplot(1,2,1)

h = pcolorRight(plotThis1,plotX,plotY); colorbar
title('log likelihood')
xlabel('mean of DFE')
ylabel('log_{10}(mutation rate)')
formatFig(gcf,gca)
set(h,'EdgeColor','none')
caxis([min(min(plotThis1)) max(max(plotThis1))])

subplot(1,2,2)
pcolorRight(log10(plotThis2),plotX,plotY); colorbar
title('log_{10}(uncertainty in log likelihood)')
caxis([min(min(log10(plotThis2))) max(max(log10(plotThis2)))])
xlabel('mean of DFE')
formatFig(gcf,gca)

%white out portions of likelihood surface some threshhold below peak
figure;
[~, thisMax] = findMaxNd(plotThis1,1);
thisMax = thisMax(1);
plotThis1( thisMax-plotThis1 > 8) = NaN;

h = pcolorRight(plotThis1,plotX,plotY); colorbar
title('log likelihood')
xlabel('mean of DFE')
ylabel('log_{10}(mutation rate)')
formatFig(gcf,gca)
set(h,'EdgeColor','none')
caxis([min(min(plotThis1)) max(max(plotThis1))])

figure;

h = pcolorRight(plotThis3,plotX,plotY); colorbar
title('log_{10}(# of simulations)')
xlabel('mean of exponential DFE')
ylabel('log_{10}(mutation rate)')
formatFig(gcf,gca)
set(h,'EdgeColor','none')

clear h aRange uRange plotThis this*

%% compute likelihood surface with patch averaging  [optional to run this cell]

%define grid spacing of DFE params in the peak region
if ~strcmp(DFE_form,'truncExp')
    deltaPar = [log10(1.2) .2*10^-2];
else
    deltaPar = [.0025 .005];
end

regBounds = zeros(length(DFEparamRange),2);
for i=1:length(DFEparamRange)
    regBounds(i,:) = [min(DFEparamRange{i}) max(DFEparamRange{i})];
end

regTicks = cell(size(DFEparamRange));
regParamRange = cell(size(DFEparamRange));
for i=1:length(DFEparamRange)
    regTicks{i} = regBounds(i,1):deltaPar(i):regBounds(i,2);
    regParamRange{i} = regTicks{i}(2:end)-deltaPar(i)/2;
end

numParams = length(DFEparamRange);
volDims = zeros(1,numParams);
for i=1:length(regTicks)
    volDims(i) = length(regTicks{i})-1;
end

patchLiks = zeros(volDims);
patchCoverage = zeros(volDims);
patchSigs = zeros(volDims);


%formulae to estimate likelihood & its uncertainty based on Wilson Score
zScore2 = 2^2; % roughly alpha=5%, see 2*(1-normcdf(zScore,0,1))
likEst = @(sawDat,simsSoFar) (sawDat + zScore2/2)./(simsSoFar + zScore2);
likSig = @(sawDat,simsSoFar) sqrt( sawDat.*(1-sawDat./simsSoFar) + zScore2/4 )./(simsSoFar + zScore2);

%linear indices of trajectory outcomes observed in the data
if trajFeature == 1
    pickTrajOutcomes = sub2ind( [size(peakFreqSdownBins,1) length(trajs)],...
            trajPeakSdown(:),...
            [1:length(trajPeakSdown)]');   
elseif trajFeature ==0
    pickTrajOutcomes = sub2ind( [2 length(trajs)],...
            trajSweepOrNot,...
            [1:length(trajs)]');       
end

for i=1:numel(patchLiks)
    
    %get the indices of parameter bounds of the box
    boxParamInds = myInd2Sub(size(patchLiks),i);
        
    %look for simulated points in the box
    pick = true(size(param2results,1),1);
    onEdge = false(size(param2results,1),1);
    for j=1:length(boxParamInds)
        pick = pick ...
            & ( param2results(:,j) > regTicks{j}(boxParamInds(j)) ...
             | abs(param2results(:,j) - regTicks{j}(boxParamInds(j))) <= 10^-numSigFigs ) ...
            & ( param2results(:,j) < regTicks{j}(boxParamInds(j)+1) ...
             | abs(param2results(:,j) - regTicks{j}(boxParamInds(j)+1))  <= 10^-numSigFigs ) ...
            ;    
    end
    
    for j=1:length(boxParamInds)
        onEdge = pick & ~(...
             param2results(:,j) > regTicks{j}(boxParamInds(j)) + 10^-numSigFigs ...
            & param2results(:,j) < regTicks{j}(boxParamInds(j)+1) - 10^-numSigFigs );
    end
        
    patchCoverage(i) = sum(pick) - (1/2)*sum(onEdge);
            
    %combine the sim results for each box
    if sum(pick) >0    
        
        weights = zeros(size(pick));
        weights(pick) = 1;
        weights(onEdge) = 1/2;
        pick = find(pick); 
        [~,~, weights] = find(weights);
        
        theseSigs = zeros(size(pick));

        theseTrajCnts = zeros([size(peakFreqSdownBins,1) length(trajs)]);
        for j=1:length(pick)
            %combine the sets of sim results within the param interval
            theseTrajCnts = theseTrajCnts + trajTypeCnts{pick(j)}.*weights(j);
        end
        theseTrajTots = sum(theseTrajCnts);
        
        if find(theseTrajTots==0,1)
            patchLiks(i) = NaN;
            patchSigs(i) = NaN;
        else 
            if trajFeature == 1
                theseTrajCnts = theseTrajCnts(pickTrajOutcomes);
            elseif trajFeature == 0
                theseTrajCnts = [theseTrajCnts(2,:); theseTrajTots(:)'-theseTrajCnts(2,:)];
                theseTrajCnts = theseTrajCnts(pickTrajOutcomes);
            end
            theseLiks = likEst(theseTrajCnts(:), theseTrajTots(:));
            theseSigs = likSig(theseTrajCnts(:), theseTrajTots(:));     
            patchLiks(i) = sum(log(theseLiks));
            patchSigs(i) = sqrt(sum( (theseSigs./theseLiks).^2));
        end
    else
        patchLiks(i) = NaN;
        patchSigs(i) = NaN;
    end    
end
clear this* these* pick* onEdge
% patchRelToPeak = patchLiks - ML;
% round(patchRelToPeak)

%%% look at patch averaged likelihood surface

plotThis1 = patchLiks(:,:);
plotThis2 = patchSigs(:,:);
plotY = regParamRange{1};
plotX = regParamRange{2};

%white out portions of likelihood surface far below peak
    thisMax = max(max(plotThis1));
     plotThis1( thisMax - plotThis1 > 5) = NaN;

h0 = figure;

subplot(1,2,1)

h1=pcolorRight(plotThis1,plotX,plotY);
colorbar
caxis([min(min(plotThis1)) max(max(plotThis1))])
formatFig(gcf,gca);
set(h1,'edgecolor','none')  

subplot(1,2,2)

h2=pcolorRight(plotThis2,plotX,plotY);
colorbar
caxis([min(min(plotThis2)) max(max(plotThis2))])
formatFig(gcf,gca);
set(h2,'edgecolor','none')
title(['uncertainty in log likelihood'])

clear plot*  h0 h1 h2
%% look at n most likely params [optional to run this cell]
n=100;
[inds topLiks] = findMaxNd(likeSurf,n);
topLikParams = zeros(size(inds,1),length(DFEparamRange));
theseSigs = zeros(size(inds,1),1);
for i=1:size(inds,1)
    for j=1:2
        topLikParams(i,j) = DFEparamRange{j}(inds(i,j));
    end
    pick = mySub2Ind(size(likeSigs),inds(i,:));
    theseSigs(i) = likeSigs(pick);
end
[topLikParams topLiks(:,1) theseSigs]

%% define peak region and check likelihoods in it have been computed densely & precisely enough [mandatory to run this cell]

%For each DFE, define peak region & "infinitessimal" grid spacing
    %DIRAC DELTA
    if strcmp(DFE_form,'dir')
             deltaPar = [log10(1.1) 1*10^-3 ];    
                peakRegBounds  = [ -6.6 -3; .02 .065];
    elseif strcmp(DFE_form,'exp')
        deltaPar = [log10(1.1) 1*10^-3 ];        
%         peakRegBounds = [ DFEparamRange{1}(peakInds(1)) + [ -2 1] ;... %min max of param 1
%                          DFEparamRange{2}(peakInds(2)) + 10^-2*[ -.5 1.5] ] ;    
        peakRegBounds = [-6 -3; 0.003 .025];
    elseif strcmp(DFE_form,'unif')
        deltaPar = [log10(1.1) 1*10^-3 ];        
        peakRegBounds  = [ -6.5 -4; .01 .05];
    elseif strcmp(DFE_form,'truncExp')
        deltaPar = [.001 .005];
            peakRegBounds = [ [0 .04]; [.06 .2]];
%        deltaPar = [log10(1.1) 10^-3 10^-5];
%          peakRegBounds = [ DFEparamRange{1}(peakInds(1)) + [ -1.5 1] ;... %min max of param 1
%                           DFEparamRange{2}(peakInds(2)) + 10^-2*[ -.5 1.25]; ...
%                           DFEparamRange{3}(peakInds(3)) + [-10^-5 10^-5] ];
    end


minLikeThresh = 11; %points with computed likelihood that's this many log likelihood units less than the ML will not be further simulated
simTop=10; %the top [simTop] most likely params will be further simulated (to beat down any params whose likelihood is high by chance)

%
paramSigFigs = ceil(log10(1./deltaPar))+1;

peakRegTicks = cell(size(DFEparamRange));
peakParamRange = cell(size(DFEparamRange));
for i=1:length(DFEparamRange)
    peakRegTicks{i} = peakRegBounds(i,1):deltaPar(i):peakRegBounds(i,2);
    peakParamRange{i} = peakRegTicks{i}(2:end)-deltaPar(i)/2;
end

numParams = length(DFEparamRange);
volDims = zeros(1,numParams);
for i=1:length(peakRegTicks)
    volDims(i) = length(peakRegTicks{i})-1;
end

%%% determine whether likelihoods in peak region have been computed
%%% densely and precisely enough
checkBoxes = zeros(volDims);
    %checkBoxes(i,j,...) = 0 => DFE param set (x,y,z, ...) needs to be
    %simulated more where (x,y,...) is approximately
    % x  =  (peakRegTicks{1}(i)+DFEparamRange{1}(i+1))/2
    % y  =  (peakRegTicks{2}(j)+DFEparamRange{2}(j+1))/2
    % etc.
    %
    
simTheseMore = cell(volDims);
boxCenterParams = cell(volDims);
    
%check the density and uncertainties of sampled points per infinitessimal vol in peak region    
boxCoverage = zeros(volDims);
    %boxCoverage(...) = number of simulated points inside vol
boxMinSig = zeros(volDims);
    %boxMinSig(...) = minimum among the uncertainties in likelihoods
    %computed for these points
boxLiks = zeros(volDims);
    %likelihood of point with minimum uncertainty in the box
boxMatchNeighbs = false(volDims);
    % = true means difference in log likelihood between this point and any 
    % of its neighbors in the peak region is less than minDelLik
minDelLik = 5; 
maxSigLike = 1/2; % %maximum allowed uncertainty
boxMinSigParams = nan(numel(boxLiks),length(DFEparamRange));
    %boxMinSigParams(i,:) = DFE params of point in box i with lowest uncertainty
boxMaxLikParams = nan(numel(boxLiks),length(DFEparamRange));
    %boxMaxLikParams(i,:) = DFE params of point in box i with max likelihood
% boxRelToPeak = boxLiks - MLE likelihood (initialized below).
    
boxSkipBCofNeighbs = zeros(volDims); 
    %if all neighbors of a param point have low likelihood, skip it.
    %boxSkipBCofNeighbs(...) = 
    %# of neighbors with likelihood < MLE - minLikeThresh
    % so if boxSkipBCofNeighbs(...) == 4, box will not be simulated
    
    %formulae to estimate likelihood & its uncertainty based on Wilson Score
    zScore2 = 2^2; % roughly alpha=5%, see 2*(1-normcdf(zScore,0,1))
    likSig = @(sawDat,simsSoFar) sqrt( sawDat.*(1-sawDat./simsSoFar) + zScore2/4 )./(simsSoFar + zScore2);
    likEst = @(sawDat,simsSoFar) (sawDat + zScore2/2)./(simsSoFar + zScore2);
%     likEst = @(sawDat,simsSoFar) sawDat./simsSoFar;
    
    %propagation of uncertainty in likelihoods to log likelihood sum
    compSig = @(theseLiks,theseSigs)  sqrt(sum( (theseSigs./theseLiks).^2)); 
    
    %linear indices of trajectory outcomes observed in the data
    if trajFeature ==1
        pickTrajOutcomes = sub2ind( [size(peakFreqSdownBins,1) length(trajs)],...
                trajPeakSdown(:),...
                [1:length(trajPeakSdown)]');     
    elseif trajFeature ==0
    pickTrajOutcomes = sub2ind( [2 length(trajs)],...
            trajSweepOrNot(:),...
            [1:length(trajs)]');     
    else error('trajFeature undefined')
    end

        
for i=1:numel(checkBoxes)
    
    %get the indices of parameter bounds of the box
    boxParamInds = myInd2Sub(size(checkBoxes),i);
        
    %look for simulated points in the box
    pick = true(size(param2results,1),1);
    onEdge = false(size(param2results,1),length(boxParamInds));
    for j=1:length(boxParamInds)
        pick = pick ...
            & ( param2results(:,j) > peakRegTicks{j}(boxParamInds(j)) ...
             | abs(param2results(:,j) - peakRegTicks{j}(boxParamInds(j))) <= 10^-numSigFigs ) ...
            & ( param2results(:,j) < peakRegTicks{j}(boxParamInds(j)+1) ...
             | abs(param2results(:,j) - peakRegTicks{j}(boxParamInds(j)+1))  <= 10^-numSigFigs ) ...
            ;
    end
        
    theseParamsInBox = param2results(pick,:);
    for j=1:length(boxParamInds)
        onEdge(:,j) = pick & ~(... %sim counts on edge between two boxes are split
             param2results(:,j) > peakRegTicks{j}(boxParamInds(j)) + 10^-numSigFigs ...
            & param2results(:,j) < peakRegTicks{j}(boxParamInds(j)+1)- 10^-numSigFigs );
    end    
    onEdge = sum(onEdge,2);
           
    %record how many points simulated in the box
    boxCoverage(i) = sum(pick) - sum((2.^(-onEdge)).*logical(onEdge));
   
    %combine the sim results within each box
    if sum(pick)>0    
        weights = zeros(size(pick));
        weights(pick) = 1;
        weights(logical(onEdge)) = 2.^(-onEdge(logical(onEdge)));
        pick = find(pick); 
        [~,~, weights] = find(weights);                
        
        theseSigs = zeros(size(pick));
        theseLiks = zeros(size(pick));
        theseTrajCnts = zeros([size(peakFreqSdownBins,1) length(trajs)]);
        
        for j=1:length(pick)
            thisSetTrajCnts = trajTypeCnts{pick(j)}*weights(j);
            thisTrajSimTot = sum(thisSetTrajCnts);
            
            %combine the sets of sim results within the param interval
            theseTrajCnts = theseTrajCnts + thisSetTrajCnts;

            %compute the likelihood and uncertainty of each set (for boxMinSig, boxMinSigParams, boxMaxLikParams)
            if trajFeature == 0
                thisSetTrajCnts = [thisSetTrajCnts(2,:); thisTrajSimTot(:)'-thisSetTrajCnts(2,:)];
            end
            thisSetSigs = likSig(thisSetTrajCnts(pickTrajOutcomes), thisTrajSimTot(:));
            thisSetLiks = likEst(thisSetTrajCnts(pickTrajOutcomes), thisTrajSimTot(:));                
            theseSigs(j) = compSig(thisSetLiks,thisSetSigs);
            theseLiks(j) = sum(log(thisSetLiks));
        end
        theseTrajTots = sum(theseTrajCnts);
        
        if find(theseTrajTots==0,1)
            boxLiks(i) = NaN;
            boxMinSig(i) = NaN;
        else
            if trajFeature == 0
                theseTrajCnts  = [theseTrajCnts(2,:); theseTrajTots(:)'-theseTrajCnts(2,:)];
            end
            theseTrajCnts = theseTrajCnts(pickTrajOutcomes);            
            boxLiks(i) = sum(log(likEst(theseTrajCnts(:), theseTrajTots(:))));
            
            [boxMinSig(i) thisSigInd] = min(theseSigs);
            [~, thisMaxLikInd] = max(theseLiks);
                 
            boxMinSigParams(i,:) = theseParamsInBox(thisSigInd,:);
            boxMaxLikParams(i,:) = theseParamsInBox(thisMaxLikInd,:);            
        end
                
    else
        boxMinSig(i) = NaN;
        boxLiks(i) = NaN;               
    end
    
end
clear these* this* pick i j
% boxCoverage
% boxMinSig

boxRelToPeak = round(boxLiks - ML)

%check differences in likelihood between nearest neighbors (only along the
%grid, not diagonal neighbors)
peakRegSize = size(boxLiks);
for i=1:numel(boxLiks)
        
   thisBoxInds = myInd2Sub(peakRegSize,i);

   %for each param axis, check the points further and back along the axis
   thisOK = true; %will be false if differences wrt any neighboring pt exceeds threshhold
   for j=1:length(thisBoxInds)
      %check behind
      if thisBoxInds(j) > 1
          thisNeighb = thisBoxInds;
          thisNeighb(j) = thisNeighb(j) - 1;          
          thisNeighb = mySub2Ind(peakRegSize,thisNeighb);
          thisOK = thisOK && abs( boxLiks(i) - boxLiks(thisNeighb)) <= minDelLik;
          
          if ML - (boxLiks(thisNeighb) + 10*boxMinSig(thisNeighb)) > minLikeThresh
              boxSkipBCofNeighbs(i) = boxSkipBCofNeighbs(i) + 1;
          end
          
      end
      
      %check forward
      if thisBoxInds(j) < peakRegSize(j)
          thisNeighb = thisBoxInds;
          thisNeighb(j) = thisNeighb(j) + 1;          
          thisNeighb = mySub2Ind(peakRegSize,thisNeighb);
          thisOK = thisOK && abs( boxLiks(i) - boxLiks(thisNeighb)) < minDelLik;
          
          if ML - (boxLiks(thisNeighb) + 10*boxMinSig(thisNeighb)) > minLikeThresh
              boxSkipBCofNeighbs(i) = boxSkipBCofNeighbs(i) + 1;
          end          
      end      
         
   end
   
   boxMatchNeighbs(i) = thisOK;   
end
 
%make list of params to simulate more
simDoMores = 0;
simDoNews = 0;
for i=1:numel(checkBoxes)
    
	thisOneDone = ...%consider sims of these params done if...
        (isfinite(boxLiks(i)) ...% likelihood is certainly far below peak
        && ML - boxLiks(i) - 10*boxMinSig(i) > minLikeThresh )... 
        || boxMinSig(i) <= maxSigLike ... %or likelihood accurately determined
        || boxSkipBCofNeighbs(i) == 2*length(DFEparamRange); %or all neighbors have low likelihood
        
    if ~thisOneDone      
        
        %if available, pick params with some simulations already     
        if isfinite(boxMinSigParams(i,1)) %lowest uncertainty
           simTheseMore{i} = [simTheseMore{i} ; boxMinSigParams(i,:)];
        end
        if isfinite(boxMaxLikParams(i,1)) %max likelihood
           simTheseMore{i} = [simTheseMore{i} ; boxMaxLikParams(i,:)];
           simDoMores = simDoMores + 1;        
        end
        %sanity check
        if ( isfinite(boxMaxLikParams(i,1)) && ~isfinite(boxMinSigParams(i,1)) ) ...
                || ( ~isfinite(boxMaxLikParams(i,1)) && isfinite(boxMinSigParams(i,1)))
            display('box has maxLikParams or minSigParams but not both? must be error')
            keyboard
        end        
        
        %if box doesn't fully contain 1 sim point, simulate params in center of box
        if boxCoverage(i) < 1  
           theseParamsToAdd = zeros(size(DFEparamRange));
           theseParamInds = myInd2Sub(peakRegSize,i);
           for j=1:length(theseParamInds)
               theseParamsToAdd(j) = ( peakRegTicks{j}(theseParamInds(j)) + peakRegTicks{j}(theseParamInds(j)+1))/2;

               %round it to the nearest sig fig
               theseParamsToAdd(j) = round(theseParamsToAdd(j)*10^paramSigFigs(j))*10^-paramSigFigs(j);

           end
           theseParamsToAdd = theseParamsToAdd(:)';

           simTheseMore{i} = [simTheseMore{i}; theseParamsToAdd];

           simDoNews = simDoNews + 1;
        end        
        
        %eliminate redundancy
        simTheseMore{i} = unique(simTheseMore{i},'rows');
    end

        checkBoxes(i) = 1;
    
end

%sim some of the most likely params more 
[inds topLiks] = findMaxNd(likeSurf,simTop);
topLikParams = zeros(size(inds,1),length(DFEparamRange));
for i=1:size(inds,1)
    for j=1:length(DFEparamRange)
        topLikParams(i,j) = DFEparamRange{j}(inds(i,j));
    end
end

for i=1:numel(checkBoxes)
    pick = true(simTop,1);
    for j=1:length(boxMaxLikParams(i,:))
        pick = pick & ...
            abs( boxMaxLikParams(i,j) - topLikParams(:,j)  ) <= 10^-numSigFigs;
    end
    if sum(pick)>0
        simTheseMore{i} = [simTheseMore{i}; boxMaxLikParams(i,:)];
        simDoMores = simDoMores + 1;
    end
end

%for reference compile a variable that stores the params at the center of
%each lattice in the peak region
for i=1:numel(boxCoverage)
    subs = myInd2Sub(volDims,i);
    boxCenterParams{i} = zeros(1,length(subs));
    for j=1:length(subs)
        boxCenterParams{i}(j) = mean([peakRegTicks{j}(subs(j)) peakRegTicks{j}(subs(j)+1)]);
    end
end

clear these* this* i j k inds
simTheseMore
simDoNews + simDoMores
%% make jobSpecs variable to launch next round of simulations
numSimSets = length(sInoc); numDims = length(DFEparamRange);

%determine size of jobSpecs var and initialize it
i = 1;
jobSpecSize = 0;
while i<= numel(simTheseMore)
    if ~isempty(simTheseMore{i})
       jobSpecSize = jobSpecSize + size(simTheseMore{i},1)*numSimSets;
        
    end    
    i = i + 1;
end
jobSpecs = zeros(1+numDims ,jobSpecSize);

jobSpecsSection = ones(1+numDims ,numSimSets); 
jobSpecsSection(1,:) = 1:numSimSets;
thisInd=1; i = 1;
while i<= numel(simTheseMore)
    if ~isempty(simTheseMore{i})
        for k=1:size(simTheseMore{i},1)
            jobSpecsSection(2:end,:) = ones(numDims,numSimSets);
            for j=1:numDims
                jobSpecsSection(j+1,:) = simTheseMore{i}(k,j);
            end

            jobSpecs(:,thisInd:thisInd+numSimSets-1) = jobSpecsSection;

            thisInd = thisInd + numSimSets;
        end
    end
    i = i+1;
end

save([DFE_form 'DFE_jobSpecs' num2str(simBatchInd+1)],'jobSpecs','DFE_form')

% Or manually specify jobSpecs variable
% simMore = .06:.002:.13;
% simThese = zeros(length(simMore),3);
% simThese(:,1) = -3.0800;
% simThese(:,2) = .008;
% simThese(:,3) = simMore(:);
% numSimSets = length(sInoc);
% 
% jobSpecs = zeros(4,size(simThese,1)*numSimSets);
% 
% thisInd = 1;
% for i=1:length(simThese)
%     pick = thisInd:thisInd+numSimSets-1;
%    jobSpecs(1,pick) = 1:numSimSets;
%    jobSpecs(2,pick) = simThese(i,1);
%    jobSpecs(3,pick) = simThese(i,2);
%    jobSpecs(4,pick) = simThese(i,3);
%     thisInd = thisInd + numSimSets;
% end

save([DFE_form 'DFE_jobSpecs' num2str(simBatchInd+1)],'jobSpecs','DFE_form')

%% WHEN ENOUGH SIMULATIONS HAVE BEEN DONE...

%ditch sim param for which not all trajectories simulated    
toss = false(size(trajTypeCnts));
for i=1:length(trajTypeCnts)
    theseTrajSimTots = sum(trajTypeCnts{i});
    if find(0==theseTrajSimTots | isnan(theseTrajSimTots) ,1)
        toss(i) = true;
    end
end
trajTypeCnts(toss) = [];
trajTypeDist(toss) = [];
trajSimTot(toss) = [];
param2results(toss,:) = [];
logLikeList(toss) = [];
logLikeSigList(toss) = [];
sum(toss)

%% record maximum likelihood and find its upper bound

dataMLE_ind = true(size(param2results,1),1);
for i=1:length(DFEparamRange)
    dataMLE_ind = dataMLE_ind & ...
        param2results(:,i) == MLE_params(i);
end
MLE_ind  = find(dataMLE_ind);

[thisMaxInds thisMax ] = findMaxNd(boxLiks,1);
thisMaxLinInd = mySub2Ind(size(boxLiks),thisMaxInds);
thisMax = boxLiks(thisMaxLinInd);
boxCenterParams{thisMaxLinInd}

% upperPeak = upBoundSurfPeak(boxLiks,peakParamRange); 
% upperPeak = upperPeak + boxMinSig(thisMaxLinInd)


% upBoundSurfPeak() sometimes fails at fminsearch(), debug later, in that case do
upperPeak = thisMax  + 3*boxMinSig(thisMaxLinInd);

clear this*
%% save likelihood surface calculations for pValFromMLE_makeWorkspace()
save(['simulation results\' ...
    [DFE_form 'DFE_compiledSims'] ],...
    'trajTypeCnts',...
    'trajTypeDist',...
    'trajSimTot',...
    'param2results',...
    'DFE_form',...
    'logLikeList',...
    'logLikeSigList',...
    'peakRegBounds',...
    'deltaPar',...
    'ML',...
    'upperPeak',...
    'MLE_ind',...
    'MLE_params');

%% to analyze a different peak, manually set ML to that peak

DFE_form = [DFE_form '2'];

n=10;
[inds topLiks] = findMaxNd(likeSurf,n);
topLikParams = zeros(size(inds,1),length(DFEparamRange));
theseSigs = zeros(size(inds,1),1);
for i=1:size(inds,1)
    for j=1:2
        topLikParams(i,j) = DFEparamRange{j}(inds(i,j));
    end
    pick = mySub2Ind(size(likeSigs),inds(i,:));
    theseSigs(i) = likeSigs(pick);
end
[topLikParams topLiks(:,1) theseSigs]

pick2ndPeak = 7;

secondPeakParams = topLikParams(pick2ndPeak,:);
secondMLind = find(param2results(:,1) == secondPeakParams(1) ...
    & param2results(:,2) == secondPeakParams(2) );
MLE_params = secondPeakParams;
MLE_ind = secondMLind;
ML = logLikeList(secondMLind);

secondPeakBox = zeros(size(MLE_params));
for i=1:length(MLE_params)
    secondPeakBox(i) = find(peakParamRange{i}>MLE_params(i),1)-1;
end
secondMaxLinInd = mySub2Ind(size(boxLiks),secondPeakBox);
boxCenterParams{secondMaxLinInd}

upperPeak = boxLiks(secondMaxLinInd)  + 3*boxMinSig(secondMaxLinInd);
%%
save(['simulation results\' ...
    [DFE_form 'DFE_compiledSims'] ],...
    'trajTypeCnts',...
    'trajTypeDist',...
    'trajSimTot',...
    'param2results',...
    'DFE_form',...
    'logLikeList',...
    'logLikeSigList',...
    'peakRegBounds',...
    'deltaPar',...
    'ML',...
    'upperPeak',...
    'MLE_ind',...
    'MLE_params');
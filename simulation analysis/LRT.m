%Likelihood ratio test: 
%is the difference in maximum likelihoods of different DFEs significant?
clear
homeDir = 'C:\Users\genya\Google Drive\BMA sim code\for distribution\';
cd(homeDir);

workspacePath = 'simulation analysis/workspaces';
addpath('simulation analysis');

%% load data and initialize functs & vars for LRT computation

load([homeDir '\lineage dynamics data\seedLinTrajs_workspaceForAnalysis.mat'])

%likelihood ratio test statistic
likeRatStat = @(ML_null,ML_alt) -(ML_null - max([ML_null ML_alt],[],2));

%peak region check
inDaBox = @(xy,boxBounds) (xy(1) >= boxBounds(1,1) && xy(1) <= boxBounds(1,2) )...
    && (xy(2) >= boxBounds(2,1) && xy(2) <= boxBounds(2,2) );

%procedural specs
trials = 2*10^2; 

trajFitClassInds = cell(length(fitClasses),1);
trajFitClassSize = zeros(length(fitClasses),1);
for i=1:length(fitClasses)
    pick = trajFitClass == fitClasses(i);
    trajFitClassInds{i} = find(pick);
    trajFitClassSize(i) = sum(pick);
end

bootTrajPeakSdown = zeros(size(trajPeakSdown));


%% load sim results and compute LRT for all pairs of DFEs

% theseDFEs = {'exp','dir','unif'};
theseDFEs = {'exp','truncExp'};
LRTpvals = zeros(length(theseDFEs));
%LRTpvals(i,j) = p-value with which null hypothesis theseDFEs{i} is
%rejected in favor of alternative hypothesis theseDFEs{j}

for i=1:length(theseDFEs)
    
    %load null hypothesis
    load([workspacePath '/' theseDFEs{i} '_simResultsForLRT.mat']);
    eval(['MLnullObs = ' theseDFEs{i} 'ML;']);
    eval(['nullMLEind = ' theseDFEs{i} 'MLE_ind;']);
    eval(['nullTrajTypeCnts = ' theseDFEs{i} 'TrajTypeCnts;']);
    eval(['nullTrajTypeDist = ' theseDFEs{i} 'TrajTypeDist;']);
    eval(['nullMLEind = ' theseDFEs{i} 'MLE_ind;']);
    eval(['nullTrajSimTot = ' theseDFEs{i} 'TrajSimTot;']);
    eval(['nullParam2Results = ' theseDFEs{i} 'Param2Results;']);
    eval(['nullPeakRegBounds = ' theseDFEs{i} '_peakRegBounds;']);
    
    %load alternative hypotheses and carry out LRT calculation
    for j=1:length(theseDFEs)
   
        if i==j
            LRTpvals(i,j) = nan;
            continue
        end
        
        load([workspacePath '/' theseDFEs{j} '_simResultsForLRT.mat']);
        
        %observed maximum likelihoods of null & alt hypotheses
        eval(['MLaltObs = ' theseDFEs{j} 'ML;']);

        if MLnullObs >= MLaltObs
            LRTpvals(i,j) = nan;
            eval(['clear ' theseDFEs{j} '*']);
            continue
        end
        
        likeRatObs = likeRatStat(MLnullObs,MLaltObs);
              
        eval(['altMLEind = ' theseDFEs{j} 'MLE_ind;']);
        eval(['altTrajTypeCnts = ' theseDFEs{j} 'TrajTypeCnts;']);
        eval(['altMLEind = ' theseDFEs{j} 'MLE_ind;']);
        eval(['altTrajSimTot = ' theseDFEs{j} 'TrajSimTot;']);
        eval(['altParam2Results = ' theseDFEs{j} 'Param2Results;']);
        eval(['altPeakRegBounds = ' theseDFEs{j} '_peakRegBounds;']);        
                
        
        %variables to store null distribution of MLs
        nullDistribLiks = zeros(trials,2);
        %nullDistribLiks(i,:) = [ML_null, ML_alt] for trial i
        peakRegCheck = zeros(trials,2);
        %peakRegCheck(i,j) = 1 => corresponding ML value in nullDistribLiks was
        %in peak region of corresponding DFE.
        
        for k=1:trials

            %generate synthetic data given the null 
            for l=1:length(bootTrajPeakSdown)
                bootTrajPeakSdown(l) = randBeans(nullTrajTypeDist{nullMLEind}(:,l),1);
            end

            %record ML value
            [nullDistribLiks(k,1), nullMLind] = findMLE(nullTrajTypeCnts,bootTrajPeakSdown,nullTrajSimTot);
            [nullDistribLiks(k,2), altMLind] = findMLE(altTrajTypeCnts,bootTrajPeakSdown,altTrajSimTot); 

            %check it's in peak region
            peakRegCheck(k,1) = inDaBox(nullParam2Results(nullMLind,:),nullPeakRegBounds);
            peakRegCheck(k,2) = inDaBox(altParam2Results(altMLind,:),altPeakRegBounds);

            if mod(k,500)==0
                display(['k = ' num2str(k)])
            end

        end
        
        if sum(sum(peakRegCheck==0))>0
%           sum(peakRegCheck)
            display('not all maximum likelihoods are in the peak region!')
        end
        
        LRTpvals(i,j) = sum(likeRatStat(nullDistribLiks(:,1),nullDistribLiks(:,2)) >= likeRatObs)/trials;
       
        eval(['clear ' theseDFEs{j} '*']);
        clear alt*
    end
    
    eval(['clear ' theseDFEs{i} '*']);
    clear null*
end

LRTpvals
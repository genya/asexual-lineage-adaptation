%% Generates workspaces for likelihood ratio test and p-value code
%takes as input sim result workspaces generated by MLEcertify.m
clear
homeDir = 'C:\Users\genya\Google Drive\BMA sim code\for distribution\';
cd(homeDir);
%%
simResultsFolder = 'simulation results';
analysisFolder = 'simulation analysis';
storeFolder = 'workspaces';

theseDFEs = {'exp','dir','unif','truncExp'}; 

for thisInd=length(theseDFEs)

    load([simResultsFolder '\' theseDFEs{thisInd} 'DFE_compiledSims.mat'])
        
    %check that if you have the right file
    if ~strcmp(DFE_form,theseDFEs{thisInd})
        keyboard
    end
    
    %save variables needed for LRT.m
    eval([theseDFEs{thisInd} 'TrajTypeCnts = trajTypeCnts;']);
    eval([theseDFEs{thisInd} 'TrajTypeDist = trajTypeDist;']);
    eval([theseDFEs{thisInd} 'TrajSimTot = trajSimTot; ']);
    eval([theseDFEs{thisInd} 'Param2Results = param2results; ']);
    eval([theseDFEs{thisInd} 'LogLikeList = logLikeList;']);
    eval([theseDFEs{thisInd} 'ML = ML;']);
    eval([theseDFEs{thisInd} 'MLE_ind = MLE_ind;']);
    eval([theseDFEs{thisInd} '_peakRegBounds = peakRegBounds;']);
    
    save([analysisFolder '\' storeFolder '\' theseDFEs{thisInd} '_simResultsForLRT'],...
        [theseDFEs{thisInd} 'ML'],...
        [theseDFEs{thisInd} 'MLE_ind'],...
        [theseDFEs{thisInd} 'TrajTypeCnts'],...
        [theseDFEs{thisInd} 'TrajTypeDist'],...
        [theseDFEs{thisInd} 'TrajSimTot'],...
        [theseDFEs{thisInd} 'Param2Results'],...
        [theseDFEs{thisInd} 'LogLikeList'],...
        [theseDFEs{thisInd} '_peakRegBounds']);
    
    %save variables needed for p-value code (pValFromMLEwrap.m)
    dataMLE_height = upperPeak;
    dataMLE_ind = MLE_ind;    
    peakRegInteriorInds = true(length(trajTypeCnts),1);
    peakRegBorderInds = true(length(trajTypeCnts),1);
    for i=1:size(peakRegBounds)
        peakRegInteriorInds  = peakRegInteriorInds ...
            & param2results(:,i) > peakRegBounds(i,1) + 2*deltaPar(i) ...
            & param2results(:,i) < peakRegBounds(i,2) - 2*deltaPar(i) ;

        peakRegBorderInds = peakRegBorderInds  ...
            & param2results(:,i) >= peakRegBounds(i,1) ...
            & param2results(:,i) <= peakRegBounds(i,2) ...
            & ~peakRegInteriorInds;    
    end
    datSize = size(trajTypeDist{1}); 
    datSize2 = [1:datSize(2)]';
    

    save([analysisFolder '\' storeFolder '\' theseDFEs{thisInd} '_simResultsForPval'],...
    'trajTypeCnts','trajSimTot','trajTypeDist',...
    'datSize', 'datSize2',...
    'dataMLE_height','dataMLE_ind',...
    'peakRegInteriorInds','peakRegBorderInds',...
    'MLE_params');
end
clear 
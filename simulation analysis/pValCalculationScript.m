clear
homeDir = 'C:\Users\genya\Google Drive\BMA sim code\for distribution\';
cd(homeDir);
simAnalysis = 'simulation analysis';
storeFolder = 'workspaces';

addpath(simAnalysis)

theseDFEs = {'exp','dir','unif'}; 

trials = 10^4;

for i = 1:length(theseDFEs)
    workspaceFilePath  = [simAnalysis '/' storeFolder '/' theseDFEs{i} '_simResultsForPval'];
    outFilePath = [simAnalysis '/pVals/'];
    
    pValFromMLEwrap(trials,theseDFEs{i},workspaceFilePath,outFilePath);%,flagFilePath);
end
    
rmpath(simAnalysis);

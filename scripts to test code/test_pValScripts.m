%function pValFromMLEwrap(trials,workspaceFile,outFile)

trials = '1000'; 
workspaceFile = 'unifDFE_simResultsForPval';
outFile = 'unifDFE_pTime';

pValFromMLEwrap(trials,workspaceFile,outFile)


%%
clear

trials = '2';
workspaceFile = 'dirDFE_simResultsForPval';
outFile = 'dirDFE_pBoot';

% profile on
synthDatPvalFromMLE_v3(trials,workspaceFile,outFile)
% profile off

% profile viewer
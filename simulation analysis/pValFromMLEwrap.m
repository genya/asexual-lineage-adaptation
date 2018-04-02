function pValFromMLEwrap(trials,DFE_form,workspaceFilePath,outFilePath,flagFilePath)
%function pValFromMLEwrap(trials,workspaceFilePath,outFilePath,flagFile)
%
%uses pValFromMLE()

load([workspaceFilePath '.mat']);

if ischar(trials)
    trials = str2double(trials);
end

[pVal, logLikeDist, passQC] = pValFromMLE(trials,...
    dataMLE_height,dataMLE_ind,...
    trajTypeDist,trajTypeCnts,trajSimTot,...
    datSize,datSize2,...
    peakRegInteriorInds,peakRegBorderInds);

if exist(outFilePath,'dir') == 0
    mkdir(outFilePath);
end

fid = fopen([outFilePath '/pVal_' DFE_form],'w');
fprintf(fid, [num2str(pVal) '\n']);
fprintf(fid, [num2str(logLikeDist(:)') '\n']);
fprintf(fid, [num2str(passQC(:)') '\n']);
fclose(fid);

%print to flagFile    
if nargin == 4
    fid = fopen([outFilePath '/pVal_' DFE_form],'a');
else
    fid = fopen([flagFilePath '/flagFile_' DFE_form],'a');
end

fprintf(fid, ['\n Fraction of p-val calcs that didn''t pass QC check 1: ' num2str(1-sum(passQC(1,:))/trials)]);
fprintf(fid, ['\n Fraction of p-val calcs that didn''t pass QC check 2: ' num2str(1-sum(passQC(2,:))/trials)]);
fprintf(fid, ['\n Fraction of p-val calcs that didn''t pass QC check 3: ' num2str(1-sum(passQC(3,:))/trials)]);
fclose(fid);




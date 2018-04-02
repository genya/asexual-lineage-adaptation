%script to generate jobSpec files & bsub submission scripts for
%trajTypeDistib()

clear
homeDir = 'C:\Users\genya\Google Drive\BMA sim code\for distribution\';
load([homeDir 'lineage dynamics data\seedLinTrajs_workspaceForSims.mat'])
%specify DFE
%     DFE_form = 'dir';
%     deltaPar = [log10(1.1) 5*10^-3]; 
%     DFEparamRange = cell(length(deltaPar),1);
%     paramBounds = [ [ -6  -3 ]  ;... %mutation rate
%                10^-2*[ .5  4 ]] ;    %mean of DFE   
 
%     DFE_form = 'exp';
%     deltaPar = [log10(2) 2.5*10^-3]; 
%     DFEparamRange = cell(length(deltaPar),1);
%     paramBounds = [ [ -7  -3 ]  ;... %mutation rate
%                10^-2*[ .5  8 ]] ;    %mean of DFE              

%     DFE_form = 'unif';
%     deltaPar = [log10(2) 10^-2]; 
%     DFEparamRange = cell(length(deltaPar),1);
%     paramBounds = [ [ -6  -3 ]  ;... %mutation rate
%                10^-2*[ .5  5 ]] ;    %mean of DFE              

%     DFE_form = 'truncExp';
    deltaPar = [.005 .005]; 
    DFEparamRange = cell(length(deltaPar),1);
    paramBounds = [ [ 0  0.05 ]  ;... %truncation from the left
                    [ 0.06  0.2 ]] ;    %truncation from the right
           %other params of distribution hard-coded to equal MLE of exponential
           % mean = 0.008
           % mutation rate = 10^-3.84

    %minimal for testing purposes
%     DFE_form = 'exp';
%     deltaPar = [log10(10) 10*10^-3]; 
%     DFEparamRange = cell(length(deltaPar),1);
%     paramBounds = [ [ -8  -3 ]  ;... %mutation rate
%                10^-2*[ 1  4 ]] ;    %mean of DFE   

for i=1:length(DFEparamRange)
    DFEparamRange{i} = paramBounds(i,1):deltaPar(i):paramBounds(i,2);
end

%specify condition(s) to exclude certain param combos
skipParam = @(params) false;
% if strcmp(DFE_form,'truncExp')
%     skipParam = @(params) params(3) <= params(2);
% else
%     skipParam = @(params) false;    
% end

%initialize job spec matrix
sInocClasses = length(sInoc);
jobSpecSize = zeros(size(DFEparamRange));
for i=1:length(DFEparamRange)
    jobSpecSize(i) = length(DFEparamRange{i});
end
jobSpecs = zeros(length(DFEparamRange)+1,prod(jobSpecSize)*sInocClasses);

%fill in job specs
thisInd = 1; thisInd2 = 1;
while thisInd <= prod(jobSpecSize) 
    
   theseInds = myInd2Sub(jobSpecSize(:)',thisInd); 
   theseParams = zeros(size(DFEparamRange));
   for i=1:length(theseInds)
       theseParams(i) = DFEparamRange{i}(theseInds(i));
   end 
   
   if sInocClasses == 1
       if skipParam(theseParams)
           jobSpecs(:,thisInd) = NaN;
       else
          jobSpecs(2:end,thisInd) = theseParams(:);
       end
   else
       if skipParam(theseParams)
           jobSpecs(:,thisInd2:(thisInd2+sInocClasses-1)) = NaN;
       else
          for j=1:length(theseParams) 
            jobSpecs(1+j,thisInd2:(thisInd2+sInocClasses-1)) = theseParams(j);
          end
          jobSpecs(1,thisInd2:(thisInd2+sInocClasses-1)) = 1:sInocClasses;
       end
   end
   
   thisInd = thisInd + 1; 
   thisInd2 = thisInd2 + sInocClasses; 
end

%clean up job spec matrix
jobSpecs(:, isnan(jobSpecs(1,:)) ) = [];

jobSpecFileName = [DFE_form 'DFE_jobSpecs0'];

save(jobSpecFileName,'jobSpecs','DFE_form')

clear this* these* i j jobSpecSize

%% generate script to run trajTypeDistrib()

simBatch = 0; %for later batches, see collect_trajTypeDistrib_outputs.m

runFile = [DFE_form 'DFE.run' num2str(simBatch)];
numJobs = size(jobSpecs,2);
serverName = 'short_serial';
errorFile = 'sim.err';
jobFolder = [DFE_form 'DFE'];

runTimeLimit = 3000; %seconds

header = ['#!/bin/sh \n'...
          '#BSUB -J "[1-' num2str(numJobs) ']" \n' ...
          '#BSUB -o /dev/null \n' ...
          '#BSUB -e ' errorFile '\n' ...
          '#BSUB -q  ' serverName '\n'];
     
runtimeTempFolder  = ['TMPDIR=/tmp/$USER${RANDOM}_${LSB_JOBINDEX} \n' ...
                      'mkdir $TMPDIR \n' ...
                      'export MCR_CACHE_ROOT=$TMPDIR \n'];
runtimeTempFolderRemove = '\n rm -rf $TMPDIR \n';
                  
programPath = ['/n/home09/efrenkel/' jobFolder '/trajTypeDistrib'];

command = [' ${LSB_JOBINDEX} ' ...
           num2str(runTimeLimit) ' '...
           DFE_form '_trajSet_${LSB_JOBINDEX} ' ...
           DFE_form '_sim' num2str(simBatch) 'Results ' ...
           jobSpecFileName];
       
mkdir(jobFolder);
fid = fopen([jobFolder '\' runFile],'w');

fprintf(fid,header);
fprintf(fid,['\n' runtimeTempFolder]);

fprintf(fid,['\n' programPath command]);

fprintf(fid,['\n' runtimeTempFolderRemove]);

fclose(fid);

%% copy code & files to job folder

filesNeededToSim = {...
    'simTraj.m',...
    'trajTypeDistrib.m',...
    'classifyTraj.m',...
    'returnGenMutFunct.m'};
for i=1:length(filesNeededToSim)
    copyfile([homeDir 'simulation code\' filesNeededToSim{i}],jobFolder);
end

copyfile([jobSpecFileName '.mat'],jobFolder);
delete([jobSpecFileName '.mat'])
copyfile([homeDir 'lineage dynamics data\' 'seedLinTrajs_workspaceForSims.mat'],jobFolder)
       


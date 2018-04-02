function trajTypeDistrib(jobInd,timeLimit,outFile,outFolder,jobSpecFile)
% function trajTypeDistrib(jobInd,timeLimit,outFile,outFolder,jobSpecFile)
%
% jobInd specifies the set of data trajectories to be simulated (via
% [jobSpecFile].mat, see below).
%
% trajTypeCnts(j,i) = # of times that simulation of data trajectory i
% yielded a trajectory of type j, 
% where j is defined by peakFreqSdownBins & classifyTraj().
%
% trajTypeCnts is written to a file 'outFolder/outFile_jobInd' with the
% results for each data trajectory written to its own line (i.e. trajTypeCnts
% is written to file with rows and columns switched).
%
% Working directory must contain: 
%       [jobSpecFile]
%       seedLinTrajs_workspaceForSims.mat
  

%hard-coded params
inocTime = 40; %when marked lineage was introduced

%timer variables
startTime = tic;
runTime = 0;
timePerSimChunk = 0;
simsPerChunk = 1;

%convert inputs to numerical data types
if ischar(jobInd)
    jobInd = str2double(jobInd);
end
if ischar(timeLimit)
    timeLimit = str2double(timeLimit);
end

%load relevant variables & parse them
load([jobSpecFile '.mat']);
simSetInd = jobSpecs(1,jobInd);
DFE_params = jobSpecs(2:end,jobInd);

load seedLinTrajs_workspaceForSims.mat
sInoc = sInoc{simSetInd};
measTimes = measTimes{simSetInd};

s = sInoc(1);
inocs = sInoc(2:end);
allMeasTimes = measTimes(:,1);
subset = logical(measTimes(:,2:end));
numDatTrajs = length(inocs);
if numDatTrajs ~= size(subset,2)
    error('number of inocula not equal to number of measurement times');
end

%obtain DFE
[genMuts delS] = returnGenMutFunct(DFE_form,DFE_params);

%seed random number generator (necessary for independence of jobs on the
%cluster)
rng(mod(cputime*now,2^30));

%initialize output
trajTypeCnts = zeros(length(peakFreqSdownBins),length(inocs));

%stopping condition variables
simChunkSize = max(floor(simsPerChunk/numDatTrajs),1)*numDatTrajs; 

%initialize temporary vars
sizTrajTypeCnts = size(trajTypeCnts);
theseTrajTypeInds = zeros(numDatTrajs,1);
datTrajInds = [1:numDatTrajs]';
theseTrajs = zeros(length(allMeasTimes),numDatTrajs);

while  runTime<timeLimit  %run until simulation time maxes out

    thisTimer = tic;

    j=0;
    while j< simChunkSize

        %simulate the set of trajs
        for i=1:numDatTrajs

            %simulate one trajectory for inoc(i)
            theseTrajs(:,i) = simTraj(s,inocs(i),allMeasTimes,genMuts,delS,inocTime);

            %take the subsets corresponding to each of the data
            %trajectory & classify them
            for k=1:numDatTrajs
                [~, theseTrajTypeInds(k)] = classifyTraj(...
                theseTrajs(:,i),allMeasTimes,[0],peakFreqSdownBins,subset(:,k));
            end
            theseTrajTypeLinInds = sub2ind(sizTrajTypeCnts,theseTrajTypeInds,datTrajInds);
            trajTypeCnts(theseTrajTypeLinInds) = trajTypeCnts(theseTrajTypeLinInds) + 1;
            
        end  

        j = j+numDatTrajs;
    end

    if timePerSimChunk == 0 || toc(thisTimer)-runTime < timePerSimChunk*3
        runTime = toc(startTime);
    else %assume that the job was suspended in the interim
        %in which case, restart the clock
        timeLimit = timeLimit - runTime;
        startTime = tic;
    end    
    
    timePerSimChunk = toc(thisTimer);

end

%print sim results
if exist(outFolder,'dir') == 0
    mkdir(outFolder);
end
fid = fopen([outFolder '/' outFile ],'w');
if fid == -1
    error('unable to open output file')
end
for i=1:numDatTrajs
    fprintf(fid,[num2str(trajTypeCnts(:,i)') '\n']);
end

%close outFile
fclose(fid);

end
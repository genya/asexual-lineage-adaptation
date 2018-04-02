%script to collect outputs of trajTypeDistribs_v8()
%Lines marked %CHECK% (there are five) are the ones that may need to be 
%updated prior to loading a new sim result set.

clear
homeDir = 'C:\Users\genya\Google Drive\BMA sim code\for distribution\';
resultsFolder = 'simulation results'; 
cd(homeDir)

%specify what simulations to load 
% simBatchInd = 23;
% DFE_form = 'unif';

% simBatchInd = 14;
% DFE_form = 'exp';

simBatchInd = 17;
DFE_form = 'dir';

% simBatchInd = 6;
% DFE_form = 'truncExp';

%location on disk of old & new workspaces
thisCompilesFolder = [homeDir '\' resultsFolder];

%name of file containing most-recent saved compilation of sim results
if simBatchInd > 0
    resultsPrev0 = [DFE_form 'DFEresults_' num2str(simBatchInd-1)];
else resultsPrev0  = '';
end

if   isempty(resultsPrev0)
    load('lineage dynamics data\seedLinTrajs_workspaceForAnalysis.mat')
else load([thisCompilesFolder '\' resultsPrev0 ]);
end
resultsPrev = resultsPrev0; clear resultsPrev0
compilesFolder = thisCompilesFolder; clear thisCompilesFolder

if ~exist('loadedFolders','var')
    loadedFolders = {};
else %display folders loaded so far
    loadedFolders(:)
end

%load data & prior sim results
resultsLatest0 = [DFE_form 'DFEresults_' num2str(simBatchInd)];

%names of folders & files specific to the sim results to-be loaded  
simFolder = [DFE_form 'DFE'];
outFolder = [DFE_form '_sim' num2str(simBatchInd) 'Results'];
jobSpecFile = [DFE_form 'DFE_jobSpecs' num2str(simBatchInd)];
clear simBatchInd %this is essential!
%%

%add this sim folder to list
loadedFolders{end+1} = [simFolder '\' outFolder];

%obtain list of DFE params & their indices
load([resultsFolder '\' simFolder '\' jobSpecFile]);

%obtain list of output files & traj set & the DFE params corresponding to each 
    %excludes empty files and file names without numerical suffixes
a = dir([resultsFolder '\' simFolder '\' outFolder]);
simFilesList = cell(length(a),1);
simFileSpecInds = zeros(size(simFilesList));
i = 1;
while i <= length(simFilesList)
    
    %does file i contain anything and have numeral at end of its name?
    r = str2double(a(i).name(end));
    if a(i).bytes > 0 && isfinite(r)
        
        %record filename
        simFilesList{i} = a(i).name;
        
        %get its numeral index
        pick = regexp(a(i).name,'_');
        simFileSpecInds(i) = str2double(a(i).name((pick(end)+1):end));
        i = i+1;
    else
        a(i) = [];
        simFilesList(i) = [];
        simFileSpecInds(i) = [];
    end
end
clear a i r pick

%intialize data structures for tabulating traj type counts
if ~exist('param2results','var') && ~exist('trajTypeCnts','var')
    DFEparamRange = cell(size(jobSpecs,1)-1,1);
    param2results = [];
    trajTypeCnts = {};
        %param2results(i,:) = DFEparams for trajTypeCnts{i}
        %the range of param j is given by unique(param2results(:,j))
end

%iterate through sim result files and add the results to trajTypeCnts
paramDim = length(DFEparamRange);
for i=1:length(simFilesList)
    
   %open sim result file
   fid = fopen([resultsFolder '\' simFolder '\' outFolder '\' simFilesList{i}],'r');
   if fid < 0
       display(['cannot open sim results file ' simFilesList{i}]);
       keyboard
   end 
   
   %determine its DFEparams and data traj indices
   thisDatTrajSet = jobSpecs(1,simFileSpecInds(i));
   theseDFEparams = jobSpecs(2:end,simFileSpecInds(i));
   
   %load the sim results
   thisNumDatTrajs = length(sim2datTraj{thisDatTrajSet});
   theseSimCnts = zeros(length(peakFreqSdownBins),thisNumDatTrajs);
   for j=1:thisNumDatTrajs
      line = fgets(fid);
      if line == -1
          display(['simulation set file ' simFilesList{i} ' does not contain results for all its trajectories'])
          break
      else
          line = str2num(line);
          theseSimCnts(:,j) = line(:);
      end
   end
   
   %close the file
   fclose(fid);
   
   %is this param combination already listed?
    listed = true(size(trajTypeCnts));
    if ~isempty(param2results)
        for j=1:paramDim
            thisASDF = theseDFEparams(j)==param2results(:,j);
            listed = listed(:) & thisASDF(:);
        end
        thisInd = find(listed);
    else thisInd = [];
    end
        
    %if yes, add it to existing results
    if ~isempty(thisInd) 
%         keyboard
        trajTypeCnts{thisInd}(:,sim2datTraj{thisDatTrajSet}) = ...
            trajTypeCnts{thisInd}(:,sim2datTraj{thisDatTrajSet}) + theseSimCnts;
    else %grow the array
        trajTypeCnts{end+1} = zeros(size(peakFreqSdownBins,1),length(trajs));
        param2results(end+1,:) = theseDFEparams(:)';
        
        trajTypeCnts{end}(:,sim2datTraj{thisDatTrajSet}) = ...
            theseSimCnts;
    end
    
end
clear fid this* i j these* line numDims listed

%update DFEparamRange
for i=1:length(DFEparamRange)
    DFEparamRange{i} = unique(param2results(:,i));
    DFEparamRange{i} = DFEparamRange{i}(:);
end
clear i

% save workspace
clear jobSpecs
clear DFEparamFile ans inFileName outFolder simFile*
resultsLatest = resultsLatest0; clear resultsLatest0;
save([compilesFolder '\' resultsLatest]);
%confidence intervals
clear
homeDir = 'C:\Users\genya\Google Drive\BMA sim code\for distribution\';
addpath('simulation analysis')
cd(homeDir);
%% generate MLE params for resampled data

simResultsFolder = 'simulation results';
analysisFolder = 'simulation analysis';
dataFolder = 'lineage dynamics data';
storeFolder = 'workspaces';

theseDFEs = {'exp','dir','unif'};

%number of trails for each DFE
trials = 10^4;

%load data
load([dataFolder '/seedLinTrajs_workspaceForAnalysis.mat'])

%initialize computation variables
randNumSamps = zeros(length(trajs),1);

%initialize output variables
MLresampParams = cell(size(theseDFEs));
MLdataParams = cell(size(theseDFEs));

for i=1:length(theseDFEs)
   
    %load sim results
    load([analysisFolder '/' storeFolder '/' theseDFEs{i} '_simResultsForLRT.mat']);
       
    %initialize output
    MLresampParams{i} = zeros(trials,2);
    eval(['MLdataParams{i} = ' theseDFEs{i} 'Param2Results(' theseDFEs{i} 'MLE_ind,:);']);
    
    for j=1:trials    
        
        %generate resampled data 
        for k=1:length(sim2datTraj)
            pick = randi(length(sim2datTraj{k}),length(sim2datTraj{k}),1);
            randNumSamps( sim2datTraj{k} ) = sim2datTraj{k}(pick);        
        end            

        %find the MLE params of the resampled data
        eval(['[~, MLresamp_ind] = findMLE( '...
            theseDFEs{i} 'TrajTypeCnts, ' ...
            'trajPeakSdown(randNumSamps(:)), '...
            theseDFEs{i} 'TrajSimTot, '...
            'randNumSamps(:));']);    
       
       eval(['MLresampParams{i}(j,:) = '...
           theseDFEs{i} 'Param2Results(' num2str(MLresamp_ind) ',:);']);
    end
       
    eval(['clear ' theseDFEs{i} '*']); 
end
clear i j k
%% plot scatter of MLE params of resampled data

%cosmetic variables
colors = 'rcbm';
symbols = 'p^sd';
markerSizes = [11 10 10 10];
lineWidths = [3 3 3 3];
xPerc = max(1/trials,.01); %plot params that come up more than xPerc % of the time
plotHandles = zeros(size(theseDFEs));
scatPlotFig = figure;

for i=1:length(theseDFEs)
    
    [theseParams, ~, theseInds] = unique(MLresampParams{i},'rows');
    toss = false(size(theseParams,1),1);
    for j = 1:size(theseParams,1)
        if sum(theseInds==j)/trials < xPerc
            toss(j) = true;
        end
    end
    theseParams(toss,:) = [];
    
    figure(scatPlotFig);

    if isempty(theseParams)
        %place holder so that legend can display
        theseParams = [nan nan];
    end
    
    plotHandles(i) = semilogy(theseParams(:,2),10.^theseParams(:,1),'o',...
                    'MarkerEdgeColor',colors(i),...
                    'MarkerFaceColor',colors(i),...
                    'MarkerSize',12);            
    hold all
   
end

for i=1:length(theseDFEs)
    
    semilogy(MLdataParams{i}(2),10.^MLdataParams{i}(1),...
        symbols(i),...
        'MarkerEdgeColor','k',...
        'MarkerSize',markerSizes(i),...
        'LineWidth',lineWidths(i));                    
end

%cosmetics
xlim([0 .07])
ylim(10.^[-7 -3])
legend(plotHandles,theseDFEs)
legend(plotHandles([1 3 2]),'exponential','uniform','\delta-function')
xlabel('mean of DFE, \rho(s)')
ylabel('mutation rate, U_b')
title([num2str(xPerc*100) '% confidence ranges for inferred DFE params'])
formatFig(gcf,gca,35)

%%
%% OVERLAY CONTOURS OF ADAPTATION RATE

load('adapLandExp.mat')

%  figure
h = gcf;
hold on
x = adapExpDFEparamRange{2};
y = adapExpDFEparamRange{1};
toss = x >.015;
x(toss) = [];
y(toss) = [];
plotThis = adapLandExp;
plotThis(toss,:) = [];
plotThis(:,toss) = [];
contour(x,...
    10.^y,...
    plotThis,...
    10^-2*[2.2],'k--','LineWidth',1.5);    
%     10^-2*[1.8 2.4],'k--','LineWidth',1.5);    
%     10^-2*[1:.5:4],'r-','LineWidth',1.5);

%%
load('adapLandDir.mat')

% figure
hold all
x = adapDirDFEparamRange{2};
y = adapDirDFEparamRange{1};
plotThis = adapLandDir;

toss = x > .04;
x(toss) = [];
plotThis(:,toss) = [];

contour(x,...
    y,...
    plotThis,...
    10^-2*[2.5],'k--','LineWidth',1.5);    


hold all
x = adapDirDFEparamRange{2};
y = adapDirDFEparamRange{1};
plotThis = adapLandDir;

toss = x < .04;
x(toss) = [];
plotThis(:,toss) = [];

contour(x,...
    y,...
    plotThis,...
    10^-2*[1.55],'k--','LineWidth',1.5);    
%     10^-2*[1:.5:4],'r-','LineWidth',1.5);

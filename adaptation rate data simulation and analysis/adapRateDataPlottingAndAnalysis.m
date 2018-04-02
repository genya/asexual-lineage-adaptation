%make figures to present adaptation rate data
clear
homeDir = 'C:\Users\genya\Google Drive\BMA sim code\for distribution';
cd(homeDir);
adapFolder = 'adaptation rate data simulation and analysis';
load([adapFolder '/' 'adaptation rate data.mat'])

addpath(adapFolder)
addpath('simulation code')
%% look at the data

postThawPropTime = 50; %generations propagated after thaw until gen20 of fitness assay
preExperimentPropTime = 40; %generations of propagation after initial clonal expansion but prior to start of experiment,
                            %should equal 'inocTime' in trajTypeDistrib.m
                            
figure;
x = [110 270 330 390 470] + ...%generations time when pops frozen down
    postThawPropTime; 

fitsOverTime = cell(length(x),1);
measTimes = zeros(size(fits));

%time point zero
ancFit = 0;
ancFitSig = std(fits(pops==0));
ancFitSem = ancFitSig/sqrt(sum(pops==0));

for i=1:16
    
   pickPops = pops == i;
   theseCols = unique(colNum(pickPops));
   plotY = zeros(size((theseCols)));
   plotErr = zeros(size((theseCols)));
   for j=1:length(theseCols)
       pickTimes = colNum == theseCols(j) & pickPops;
       plotY(j) = mean(fits(pickTimes));
       plotErr(j) = std(fits(pickTimes));
       measTimes(pickTimes) = x(j);
       fitsOverTime{j} = [fitsOverTime{j}; fits(pickTimes)];
   end
      
    errorbar([0; x(:)],[ancFit; plotY(:)],[ancFitSig; plotErr(:)],'.-','LineWidth',2,'Color',[.9 .9 .9]);
    hold all
    
end
 ylim([-.02 .17])

xlabel('generations')
ylabel('fitness')
% xlim([0 max(measTimes)])
formatFig(gcf,gca)

fitMean = zeros(size(x));
fitStd = zeros(size(x));
for i=1:length(x)
    fitMean(i) = mean(fitsOverTime{i});
    fitStd(i) = std(fitsOverTime{i});
end
hold all
errorbar([0; x(:)],[ancFit; fitMean(:)],[ancFitSig; fitStd(:)],'.-','LineWidth',2,'Color','k')


% compile data into sensible variables

popTrajMean = zeros(length(x),16); 
popTrajStd = zeros(length(x),16);
popTrajSem = zeros(length(x),16);
popTrajRawDat = cell(16,1);

for i=1:16
    
   popTrajRawDat{i} = nan(4,length(x));
    
   for j=1:length(x)
       pick = measTimes == x(j) & pops == i;
       
       popTrajRawDat{i}(1:sum(pick),j) = fits(pick);
       
       popTrajMean(j,i) = mean(fits(pick));
       popTrajStd(j,i) = std(fits(pick));       
       popTrajSem(j,i) = popTrajStd(j,i)/sum(pick);       
   end
end
popTrajVar = popTrajStd.^2;

popTrajAdapRate = zeros(length(x)-1,16);
popTrajAdapStd = zeros(length(x)-1,16);
for i=1:length(x)-1
   popTrajAdapRate(i,:) =  10^2*(popTrajMean(i+1,:)-popTrajMean(i,:))/(x(i+1)-x(i));
   popTrajAdapStd(i,:)  =  10^2*sqrt(popTrajVar(i+1,:)+popTrajVar(i,:))/(x(i+1)-x(i));
end

%progate error to meanAdapRate
popFitTotMean = mean(popTrajMean,2);
popFitTotStd = sqrt(  sum(popTrajStd.^2,2));

%add zero time point
popFitTotMean = [ancFit; popFitTotMean];
popFitTotStd = [ancFitSig; popFitTotStd];
popFitTotSem = popFitTotStd/sqrt(16);
x = [0 x];
popTrajStd = [ ones(1,size(popTrajStd,2))*ancFitSig; popTrajStd ];

clear i j pick* plot* this* these*

%% simulate the data

%sim inputs & ouputs (MLE params for each DFE)
DFE_forms = {'dir','exp','unif'};
DFE_params = {...
    [-6.29 0.0585],...
    [-4.00 .0085],...
    [-5.69 0.032]};

measTimes = x;
measTimesSmooth = min(x):10:max(x);
measPoints = zeros(size(measTimes));
for i=1:length(x)
    [~, measPoints(i)] = min(abs(measTimesSmooth-measTimes(i)));
end

batchSize = 16;
trials = 250;

predAdapRate = cell(size(DFE_forms));
predAdapRateSem = cell(size(DFE_forms));

predAdapRateSmooth = cell(size(DFE_forms));
predAdapRateSemSmooth = cell(size(DFE_forms));

%simulate
for k=1:length(DFE_forms)

    [genMuts delS]  = returnGenMutFunct(DFE_forms{k},DFE_params{k});

    simTrajs = zeros(length(measTimes),trials);
    batchTrajsSmooth = zeros(length(measTimesSmooth),batchSize);
    simTrajsSmooth = zeros(length(measTimesSmooth),trials);

    for i=1:trials
        for j=1:batchSize
            batchTrajsSmooth(:,j) = adapRate(measTimesSmooth + preExperimentPropTime,genMuts,delS);        
        end

        %sample time points measured
        batchTrajs = batchTrajsSmooth(measPoints,:);
        
        %add meas noise 
%         batchTrajs = batchTrajs  + randn(size(popTrajStd)).*popTrajStd;
        batchTrajs = batchTrajs  + randn(size(popTrajStd)).*popTrajStd;

        simTrajsSmooth(:,i) = mean(batchTrajsSmooth,2);
        simTrajs(:,i) = mean(batchTrajs,2);
    end

    predAdapRate{k} = mean(simTrajs,2);
    predAdapRateSem{k} = std(simTrajs,0,2); %standard deviation of null distribution
    predAdapRateSmooth{k} = mean(simTrajsSmooth,2);
    predAdapRateSemSmooth{k} = std(simTrajsSmooth,0,2); %standard deviation of null distribution
end
clear k genMuts delS i j

%% plot

%cosmetic vars
h1 = figure;
xlabelText = 'Generations';
ylabelText = 'Fitness (%)';
fs = 18;
curveHandles = zeros(2+length(DFE_forms),1);
curveHandlesSmooth = zeros(length(DFE_forms),1);

%pop avg
curveHandles(1) = errorbar(x,popFitTotMean*100,popFitTotSem*100,'k.-','LineWidth',4,'MarkerSize',20);

%raw dat
% figure(h1);
% for i=1:size(popTrajMean,2)    
%     hold all 
% %     this = plot(x(:),...
% %         100*[ancFit; popTrajMean(:,i)],'.-','LineWidth',2,'Color',[.9 .9 .9]);
%     this = errorbar(x(:),...
%         100*[ancFit; popTrajMean(:,i)],100*[ancFitSem; popTrajSem(:,i)],'.-','LineWidth',2,'Color',[.9 .9 .9]);
%     uistack(this,'bottom')
% end
% curveHandles(2) = this;

%sim avg
colors = 'crb';
lineStyle = {'--','--','--'};
for i=1:length(DFE_forms)
    hold all
    curveHandles(2+i) = errorbar(x(:), predAdapRate{i}(:)*100, predAdapRateSem{i}(:)*100,...
        [colors(i) '.'],'LineWidth',3);
    
    curveHandlesSmooth(i) = plot(measTimesSmooth, predAdapRateSmooth{i}(:)*100,...
        [colors(i) lineStyle{i}],'LineWidth',3);
    if strcmp(DFE_forms,'exp')
        uistack(curveHandlesSmooth(i),'top')
        uistack(curveHandles(2+i),'top')
    end
end


%cosmetics
% formatFig(gcf,gca)
set(gca,'FontSize',fs)
xlabel(xlabelText,'FontSize',fs)
ylabel(ylabelText,'FontSize',fs)
set(gca,'YTick',[0:2:10])
xlim([0 max(x)+10])
ylim([-.5 10])

%reorder legend
% thisLegCurves = [curveHandles(1); curveHandles(4); curveHandles(3); curveHandles(5)];%; curveHandles(2)];
thisLegCurves = [curveHandles(1); ...
    curveHandlesSmooth(2); curveHandlesSmooth(3); curveHandlesSmooth(1)];
a=legend(thisLegCurves,...
    'average population fitness',...
    'exponential DFE prediction',...
    'uniform DFE prediciton',...
    '\delta-DFE prediction',... %,'individual population fitnesses',...
    'Location','NorthWest');
set(a,'FontSize',17)

%% plot SI figure
%cosmetic vars
h1 = figure;
xlabelText = 'Generations';
ylabelText = 'Fitness (%)';
fs = 18;
curveHandles = zeros(2,1);

%pop avg
curveHandles(1) = errorbar(x,popFitTotMean*100,popFitTotSem*100,'k.--','LineWidth',2,'MarkerSize',20);
% curveHandles(1) = plot(x,popFitTotMean*100,'k--','LineWidth',2,'MarkerSize',10);

%raw dat
figure(h1);
for i=1:size(popTrajMean,2)    
    hold all 
%     this = plot(x(:),...
%         100*[ancFit; popTrajMean(:,i)],'.-','LineWidth',2,'Color',[.9 .9 .9]);
    this = errorbar(x(:),...
        100*[ancFit; popTrajMean(:,i)],100*[ancFitSem; popTrajSem(:,i)],'.-','LineWidth',2,'Color',[.7 .7 .7]);
    uistack(this,'bottom')
end
curveHandles(2) = this;

%cosmetics
% formatFig(gcf,gca)
set(gca,'FontSize',fs)
xlabel(xlabelText,'FontSize',fs)
ylabel(ylabelText,'FontSize',fs)
set(gca,'YTick',[0:2:16])
xlim([-10 max(x)+10])
ylim([-.5 16])

%reorder legend
% a=legend(curveHandles,'average population fitness','individual population fitnesses',...
%     'Location','NorthWest');
% set(a,'FontSize',17)


%% is the adaptation rate before and after gen 300 significantly different?
% look at the distribution of this difference for bootstrapped data.

trials = 10^4;
randDraws = randi(size(popTrajMean,2),trials,size(popTrajMean,2));

timeIntervals = zeros(4,16);
for i=2:length(x)-1
    timeIntervals(i-1,:) = (x(i+1)-x(i))/100;
end

adapRateDist = zeros(trials,2); 
  %adapRateDist(i,:) = [mean rate of adap for gen <= 300; for gen >= 300]
  
k = 1;  
for i=1:trials
%    popTrajAdapDraw = popTrajAdapRate(:,randDraws(i,:)) + randn(size(popTrajAdapStd)).*popTrajAdapStd;
                            %change this to an empirical bootstrap, draw
                            %four measurements with raplcement and take the mean 
                          
                            
   randDat = randi(4,5,length(randDraws(i,:)));                         
   for j=1:length(randDraws(i,:))                        
        pick = sub2ind(size(popTrajRawDat{j}),randDat(:,j),[1:5]');
        randDat(:,j) = popTrajRawDat{j}(pick);        
   end
   popTrajAdapDraw = [randDat(2:end,:)-randDat(1:end-1,:) ]./timeIntervals;                          
                
   adapRateDist(i,:) = [mean(mean(popTrajAdapDraw(1:2,:),2)) mean(mean(popTrajAdapDraw(3:4,:),2))];
end
% adapRateDist 

[n1 x1] = hist(adapRateDist(:,1),20);
[n2 x2] = hist(adapRateDist(:,2),20);
[n3 x3] = hist(adapRateDist(:,1)-adapRateDist(:,2),20);
[n4 x4] = hist(adapRateDist(:,2)./adapRateDist(:,1),20);
figure
subplot(2,1,1)
plot(x1,n1,'.-',x2,n2,'.-',x3,n3,'.-')


xlabel('adaptation rate v (% per 100 generations)')
ylabel('times drawn by bootstrap')
legend('v before gen 300', 'v after gen 300','v_{<300}-v_{>300}','Location','NorthOutside','Orientation','horizontal')
formatFig(gcf,gca)

subplot(2,1,2)
plot(x4,n4,'.-')
xlabel('v_{>300}/v_{<300}')
ylabel('times drawn by bootstrap')
formatFig(gcf,gca)

deltaVobs = [mean(mean(popTrajAdapRate(1:2,:))) mean(mean(popTrajAdapRate(3:4,:)))];
deltaVobs = deltaVobs(2)-deltaVobs(1);

deltaVnull = deltaVobs  - (adapRateDist(:,2)-adapRateDist(:,1));
pVal = sum(abs(deltaVnull)>abs(deltaVobs))/trials

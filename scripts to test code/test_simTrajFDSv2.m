% script to test simTraj()
%% look at trajectories

%DFE specs
    %for single param DFEs
    DFE_form = 'exp';
%     DFE_params = [-4 .0085];
%       DFE_params = [-4+log10(2) .0085];
      DFE_params = [-4 .0085];
%     DFE_params = [-4 .005];
    
%     DFE_form = 'truncExp';
%     DFE_params = [-4 .0085 .02 .073];

%     DFE_form = 'dir';
%     DFE_params = [-4 .02];
   
%experiment params
measTimes = [0:10:10^4];


trials = 12;
initNomoEqu=linspace(.1,.9,trials);

trajs = zeros(length(measTimes),trials);

[genMuts delS] = returnGenMutFunct(DFE_form,DFE_params);

for i=1:trials
 	trajs(:,i) = simTrajFDSv2(measTimes,genMuts,delS,initNomoEqu(i));
%     trajs(:,i) = simTrajNoFDS(measTimes,genMuts,delS,initNomoEqu);
end
% xlim([0 5000])
mean(trajs(end,:))


figure
pickTimes = 1:3:length(measTimes);
plot(measTimes(pickTimes),trajs(pickTimes,:),'k-','LineWidth',3)
ylim([0 1])
xlim([0 5000])
formatFig(gcf,gca,19)

% for i=1:12
%     subplot(3,4,i)
%     % detGrow = (inoc/(1-inoc))*exp(s*(measTimes));
%     % detGrow = detGrow./(1+detGrow);
%     pickTimes = 1:4:length(measTimes);
%     pickTrajs = [floor((i-1)*trials/12) + 1 : floor(i*trials/12)];
%     plot(measTimes(pickTimes),trajs(pickTimes,pickTrajs),'k-','LineWidth',3)
%     % hold all
%     % plot(measTimes,detGrow,'k--','LineWidth',3)
%     %xlabel('generations')
%     %ylabel('frequency seeded lineage')
%     title(num2str(i))
%     ylim([0 1])
%     %xlim([400 1.4*10^3])
%     % xlim([0 5*10^3])
%     formatFig(gcf,gca,19)
% end
% title('\DeltaS = 0.25')


%%

DFE_form = 'exp';
DFE_params = [-4 .0085];
[genMuts delS] = returnGenMutFunct(DFE_form,DFE_params);

startFreqs = .1:.2:.9;%0:.025:1;
%startFreqs = [10.^[-3:.1:-1] .15:.05:.95 .99 1];
pFixVsFreq = zeros(size(startFreqs));
time2FixVsFreq = zeros(length(startFreqs),2);
time2FixVsFreqStd = zeros(length(startFreqs),2);

pFixVsFreqNoFDS = zeros(size(startFreqs));
time2FixVsFreqNoFDS = zeros(length(startFreqs),2);
time2FixVsFreqStdNoFDS = zeros(length(startFreqs),2);

trials = 50;%2*10^3;
measTimes = [0:10:10^5];
%%
for j=1:length(startFreqs)
    
    
    %with FDS
    theseTimes2FixA = [];
    theseTimes2FixB = [];
    
    for i=1:trials
       thisTraj = simTrajFDSv2(measTimes,genMuts,delS,startFreqs(j));
%         thisTraj = simTrajFDSv3(measTimes,genMuts,delS,startFreqs(j));
        
        if thisTraj(end) == 1
            pFixVsFreq(j) = pFixVsFreq(j) + 1;            
            theseTimes2FixA = [theseTimes2FixA; measTimes(find(thisTraj==1,1))];
        elseif thisTraj(end) == 0
            theseTimes2FixB = [theseTimes2FixB; measTimes(find(thisTraj==0,1))];
        else
            display(['traj unresolved, startFreq = ' num2str(startFreqs(j))])
        end    
    display(['with FDS trial ' num2str(i)])
    end     
    
    time2FixVsFreq(j,1) = mean(theseTimes2FixA);
    time2FixVsFreqStd(j,1) = std(theseTimes2FixA);
    time2FixVsFreq(j,2) = mean(theseTimes2FixB);
    time2FixVsFreqStd(j,2) = std(theseTimes2FixB);    
    
    %no FDS
%     theseTimes2FixA = [];
%     theseTimes2FixB = [];
%     
%     for i=1:trials
%        thisTraj = simTrajNoFDS(measTimes,genMuts,delS,startFreqs(j));
%        
%         if thisTraj(end) == 1
%             pFixVsFreqNoFDS(j) = pFixVsFreqNoFDS(j) + 1;            
%             theseTimes2FixA = [theseTimes2FixA; measTimes(find(thisTraj==1,1))];
%         elseif thisTraj(end) == 0
%             theseTimes2FixB = [theseTimes2FixB; measTimes(find(thisTraj==0,1))];
%         else
%             display(['non-FDS traj unresolved, startFreq = ' num2str(startFreqs(j))])
%         end    
%     display(['no FDS trial ' num2str(i)])
%     end     
%     
%     time2FixVsFreqNoFDS(j,1) = mean(theseTimes2FixA);
%     time2FixVsFreqStdNoFDS(j,1) = std(theseTimes2FixA);
%     time2FixVsFreqNoFDS(j,2) = mean(theseTimes2FixB);
%     time2FixVsFreqStdNoFDS(j,2) = std(theseTimes2FixB);        
    
    j
end

pFixVsFreq = pFixVsFreq/trials;
pFixVsFreqNoFDS = pFixVsFreqNoFDS/trials;

%% plot

fontSize = 19;
lineWidth = 3;
markerSize = 25;

figure
% subplot(1,2,1)
plot(startFreqs,pFixVsFreq,'r.-','LineWidth',lineWidth,'MarkerSize',markerSize)
hold on 
% plot(startFreqs,pFixVsFreqNoFDS,'k.','LineWidth',lineWidth,'MarkerSize',20)
xlim([0 1])
plot([0 1], [0 1],'k--','LineWidth',lineWidth)
formatFig(gcf,gca,fontSize)
ylim([-.01 1.05])

figure
% subplot(1,2,2)
plot(startFreqs,time2FixVsFreq,'r.-','LineWidth',lineWidth,'MarkerSize',markerSize)
%errorbar(startFreqs,time2FixVsFreq,time2FixVsFreqStd/sqrt(trials-1),'b.-')
hold on
plot(startFreqs,time2FixVsFreqNoFDS,'k.-','LineWidth',lineWidth,'MarkerSize',markerSize)
%errorbar(startFreqs,time2FixVsFreqNoFDS,time2FixVsFreqStdNoFDS/sqrt(trials-1),'k.-')
xlim([0 1])
formatFig(gcf,gca,fontSize)
vline(0.5)

figure
% subplot(1,2,2)
plot(startFreqs,time2FixVsFreq(:,1),'r.-','LineWidth',lineWidth,'MarkerSize',markerSize)
%errorbar(startFreqs,time2FixVsFreq,time2FixVsFreqStd/sqrt(trials-1),'b.-')
hold on
plot(startFreqs,time2FixVsFreqNoFDS(:,1),'k.-','LineWidth',lineWidth,'MarkerSize',markerSize)
%errorbar(startFreqs,time2FixVsFreqNoFDS,time2FixVsFreqStdNoFDS/sqrt(trials-1),'k.-')
xlim([0 1])
formatFig(gcf,gca,fontSize)


figure
% subplot(1,2,2)
plotThis = time2FixVsFreq(:,1).*pFixVsFreq(:) + time2FixVsFreq(:,2).*(1-pFixVsFreq(:));
plot(startFreqs,plotThis ,'r.-','LineWidth',lineWidth,'MarkerSize',markerSize)
hold on
plotThis = time2FixVsFreqNoFDS(:,1).*pFixVsFreqNoFDS(:) + time2FixVsFreqNoFDS(:,2).*(1-pFixVsFreqNoFDS(:));
plot(startFreqs,plotThis,'k.-','LineWidth',lineWidth,'MarkerSize',markerSize)
xlim([0 1])
ylim([0 2600])
formatFig(gcf,gca,fontSize)

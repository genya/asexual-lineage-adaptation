%test classifyTraj() by manually defining seeded lineage trajectories and
%seeing that they are properly classified

load('seedLinTrajs_workspaceForSims.mat')
clear measTimes sInoc

%manually define seeded lineage trajectory
traj = [.01 .05 .1 .2 .1 .01 .5 .01 0];
measTimes = [0 100*(1:(length(traj)-1))];
indOfFirstPeak = 4;

%fitness after peak
trajRat = traj./(1-traj);
sTrue = log(trajRat(indOfFirstPeak)/trajRat(indOfFirstPeak+1))/(measTimes(indOfFirstPeak+1)-measTimes(indOfFirstPeak));

pick = true(size(traj));
xtremaBins = [0 1 2 3 4];

[xtremaBinned peakFreqSdownBin] = classifyTraj(traj,measTimes,xtremaBins, peakFreqSdownBins,pick);


%plot
figure
plot(measTimes,traj,'.-')
title(['number of extrema = ' num2str(xtremaBins(xtremaBinned))])
ylim([-.1 1.1])

%show peak sector
peakSectors = unique(peakFreqSdownBins(:,1));
peakSectors = peakSectors(3)-peakSectors(2);
peakSectors = 0:peakSectors:1;

hold all
peakFreqSdownBins(peakFreqSdownBin,1)
pick = [find(peakSectors <  peakFreqSdownBins(peakFreqSdownBin,1),1,'last'),...
    find(peakSectors >  peakFreqSdownBins(peakFreqSdownBin,1),1,'first')];
hline(peakSectors(pick(1)));
hline(peakSectors(pick(2)));

display('binned (f_{peak},s_{down}) and actual (f_{peak},s_{down})')
[peakFreqSdownBins(peakFreqSdownBin,:);
    traj(indOfFirstPeak) sTrue]

peakFreqSdownBins

%% test classifyTraj() by inspecting its classification of data trajectories


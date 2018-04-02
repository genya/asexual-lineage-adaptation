% script to test simTraj()
%% look at trajectories

%DFE specs
    %for single param DFEs
%     DFE_form = 'exp';
%     DFE_params = [-4 .0085];

%     DFE_form = 'dir';
%     DFE_params = [-4 .02];
   
    DFE_form = 'truncExp';
    DFE_params = [-4 .0085 .02 .073];

%experiment params
s = 0;
measTimes = [0:10:10^3];
inoc = 0.5;

inocTime = 10;

trials = 10;

trajs = zeros(length(measTimes),trials);

[genMuts delS] = returnGenMutFunct(DFE_form,DFE_params);

for i=1:trials
	trajs(:,i) = simTraj(s,inoc,measTimes,genMuts,delS,inocTime);
end

%% 
figure
detGrow = (inoc/(1-inoc))*exp(s*(measTimes));
detGrow = detGrow./(1+detGrow);

plot(measTimes,trajs,'.-')
hold all
plot(measTimes,detGrow,'k--','LineWidth',3)
xlabel('generations')
ylabel('frequency seeded lineage')


%% does simTraj() correctly simulate genetic drift?

%number of simulations per (inoc,s_seed)
trials = 100;

%blank DFE
[genMuts delS] = returnGenMutFunct('dir',[-10 0]);

%establishment probability
Nb = 10^4;
Ne = Nb;
kimura = @(s,inoc) (1 - exp(-2*10*s*Ne*inoc) )...
    ./(1 - exp(-2*10*s*Ne));  %standard Kimura formula with s-->10*s because of the 10 generations between bottlenecks

%inoc params
s_seed = 10.^linspace(log10(0.001),log10(.075),6);
inoc = 10.^linspace(log10(1/Nb),log10(.01),6);
inocParams = zeros(length(s_seed)*length(inoc),2);
for i=1:length(s_seed)
    for j=1:length(inoc)
        inocParams( (i-1)*length(inoc) + j,:) = [s_seed(i) inoc(j)];
    end
end

%analytical prediction
%toss params for which analytical prediction of drift escape >99% or < 1%
driftPred = zeros(size(inocParams,1),1);
i = 1;
while i<=length(driftPred)
    thisProb = kimura(inocParams(i,1),inocParams(i,2));
    if thisProb > .99 || thisProb < .01
        inocParams(i,:) = [];
        driftPred(i) = [];
    else
        driftPred(i) = thisProb;
        i = i+1;
    end
end

%simulated drift escape
driftSim = zeros(size(driftPred));
measTimes = [0:50:5*10^3];
% profile on
for i=1:size(inocParams,1)
    trajs = zeros(length(measTimes),trials);
    for k=1:trials
       trajs(:,k) = simTraj(...
           inocParams(i,1),...
           inocParams(i,2),...
           measTimes,genMuts,delS,0);
    end
    driftSim( i) =...
        sum(trajs(end,:)>=.1)/trials;
end
figure
errorbar(driftPred,driftSim,sqrt(driftPred/trials),'.')
hold on
plot(driftPred,driftPred,'k--')
hold off
xlabel('\pi(s,inoc) = kimura formula')
ylabel('\pi(s,inoc) from simulations')
title([num2str(trials) ' simulations per point'])

figure
sProdInoc = inocParams(:,1).*inocParams(:,2)*Nb;
[~, inds] = sort(sProdInoc);
plot(log(sProdInoc(inds)),driftSim(inds),'.',...
    log(sProdInoc(inds)),driftPred(inds),'k--')
hold on
errorbar(log(sProdInoc(inds)),driftSim(inds),sqrt(driftPred(inds)/trials),'.')
xlabel('N_b*log(inoc*s)')
ylabel('probability to escape drift')
legend('simulations','kimura \pi(s,inoc)','Location','NorthWest')
hold off

%Define DFE
%  DFE_form = 'exp';
% DFE_params = [-3.14 .0081];

% DFE_form = 'dir';
% DFE_params = [-4.03 .0381];

DFE_form = 'unif';
DFE_params = [-4 .0495];


%fitness & inoc of seeded lineage
s = .04;
inoc = 0.0143;

%duration and number of sims
measTimes = [0:20:600];
trials = 10;
Nb = 10^4;
inocTime = 40;

trajs = zeros(length(measTimes),trials);

[genMuts delS] = returnGenMutFunct(DFE_form,DFE_params);

for i=1:trials
   trajs(:,i) = simTraj(s,inoc,measTimes,genMuts,delS,inocTime );
end

figure
plot(measTimes,trajs,'.-')

determGrowth = (inoc/(1-inoc))*exp(s*measTimes);
determGrowth = determGrowth./(1+determGrowth);

hold on
plot(measTimes, determGrowth,'k--','LineWidth',3)
%

%sweep or not
% thresh = .99;
% pick = find(determGrowth>thresh,1);
% if isempty(measTimes(pick))
%     pick = length(measTimes);
% end
% pSweep = sum(trajs(pick,:)>.95)/trials

% pSweep = 0;
% for i=1:trials
%     [maxtab mintab] = peakdet(trajs(:,i),.01);
%     if length(maxtab) == 0
%         pSweep = pSweep + 1;
%     end
% end
% pSweep = pSweep/trials

pSweep = sum(trajs(end,:)>.5)/trials

title([DFE_form ' DFE, s_{seed} = ' num2str(s) ', inoc = ' num2str(inoc*100) '% , N_b = 10^{' num2str(log10(Nb)) '}, p_{sweep} \approx ' num2str(pSweep)])
%%


%% generate synthetic data figure
clear
addpath('C:\Users\genya\Google Drive\BMA sim code')
load('C:\Users\genya\Google Drive\BMA sim code\code for analyzing outputs\BMAdata_workspaceForSimAnalysis.mat')
%%
trials = 0;
s = fitClasses;
thresh = .05;

simNotPickRandTrajs = false;
figure
%MLE stuff
%     DFE_form = 'exp';
%     DFE_params = [-3.3 .0085];
    DFE_form = 'exp';
    DFE_params = [-3.8 .008];


%     DFE_form = 'dir';
%     DFE_params = [-4.5 .032];

%     DFE_form = 'dir';
%     DFE_params = [-6.25 .055];

%     DFE_form = 'exp';
%     DFE_params = [-3.83+.5 .008];
    
%     DFE_form = 'truncExp';
%     DFE_params = [-4 .01 .073];
    
    
% DFE_form = 'piecewise';
% DFE_params = [-10.0000   -6.8650   -6.3100   -6.8525   -7.5450];

if simNotPickRandTrajs
    [genMuts delS] = returnGenMutFunct2(DFE_form,DFE_params);    
end
figure
%cosmetics
cmap = colormap; 
cmap = get(gcf,'Colormap');
close(gcf);
colorScale = linspace(...
        -3,...
        -2,...
        size(cmap,1));



%for each fitness class
   figure
for i=1:length(s) 
    
  
  pickThese = find(trajFitClass  == fitClasses(i));
  if trials > 0
      pickThese = pickThese(randi(length(pickThese),[trials 1]));
  end

  %time range
  thisXaxis = [0 0];
   
  subplot(2,2,i)
% subplot(3,3,5)
  
  for j=1:length(pickThese)
      
    theseTimes =  trajs{pickThese(j)}(1,:)*10; 
    thisInocYR =  trajs{pickThese(j)}(2,1);    
    thisInocFreq = inocFreqBestFit(pickThese(j));
    
    
    if simNotPickRandTrajs

        thisTraj = simTraj13withSGV2(s(i),...
            thisInocFreq, ...
            theseTimes,...
            genMuts,delS,...
            30);        
        
%         thisTraj = simTraj13(s(i),...
%             thisInocFreq, ...
%             theseTimes,...
%             genMuts,delS);
    else
     thisTraj = trajs{pickThese(j)}(2,:)./(trajs{pickThese(j)}(2,:)+1); 
    end

    %identify when traj crossed thresh
    timeInd2 = find(thisTraj>=thresh,1,'first');
    timeInd1 = find(thisTraj(1:timeInd2)<=thresh ,1,'last');
    timeInds  =  unique([timeInd1 timeInd2 ]); clear timeInd1 timeInd2
    
    %interpolate to find moment when thresh crossed
    theseY = thisTraj(timeInds);
    theseX = theseTimes(timeInds);
    if length(timeInds)==2
        
        %lineare interpolation
%         slope =  ( theseY(2) - theseY(1))/(theseX(2)-theseX(1));
%         %solve: theseY(1) + deltaT*slope = thresh
%         deltaT = (thresh - theseY(1))/slope + theseX(1);
            %same thing but less transparent:        
            %         deltaT = ( theseX(1)*(theseY(2)-thresh) + theseX(2)*(thresh - theseY(1)) )...
            %             /( theseY(2) - theseY(1));
        
        %exponential interpolation
        theseRats = theseY./(1-theseY);
        threshRat = thresh/(1-thresh);
        thisFit = log( theseRats(2)/theseRats(1))/(theseX(2)-theseX(1));
%         thisFit = s(i);
        %solve for detalT: threshRat = theseRats(1)*exp(thisFit*(deltaT-theseX(1))
        deltaT = log(threshRat/theseRats(1))/thisFit + theseX(1);
        
    elseif length(timeInds) == 1
        if thisTraj(1) > thresh
            deltaT = 0;
        else
            keyboard
        end
    else
        
        %find a point within 1% of thresh
        pick = find(abs(thisTraj-thresh)<.01,1,'first');
        if ~isempty(pick)
            deltaT = theseTimes(pick);
        else keyboard
        end
%         thisT = (1/s(i))*log((thresh/(1-thresh))/inocFreqBestFit(pickThese(j)));
    end

    theseTimes = theseTimes-deltaT;    
    thisXaxis = [min(thisXaxis(1),theseTimes(1)) max(thisXaxis(2),theseTimes(end))];
    
    
    thisColor = logical(hist(log10(thisInocYR),colorScale));
    
%     %interpolate each traj
    thisTrajInterpolTimes = theseTimes(1):10:theseTimes(end);
    thisTrajInterpol = zeros(size(thisTrajInterpolTimes ));
    for k=1:length(thisTrajInterpol)
        
        %find surrounding time points
        timeInd2 = find(theseTimes>=thisTrajInterpolTimes(k),1,'first');
        timeInd1 = find(theseTimes(1:timeInd2)<=thisTrajInterpolTimes(k),1,'last');
        timeInds  =  unique([timeInd1 timeInd2 ]); clear timeInd1 timeInd2
        
        if length(timeInds) == 2
            theseY = thisTraj(timeInds);
            theseX = theseTimes(timeInds);        
            theseRats = theseY./(1-theseY);
            thisFit = log( theseRats(2)/theseRats(1))/(theseX(2)-theseX(1));        
            
            thisRat = theseRats(1)*exp(thisFit*(thisTrajInterpolTimes(k) - theseX(1)));
            thisTrajInterpol(k) = thisRat./(1+thisRat);
            
        elseif length(timeInds) == 1
            thisTrajInterpol(k) = thisTraj(timeInds);
        else 
            keyboard
        end       
    end
    
    plot(thisTrajInterpolTimes,thisTrajInterpol,'-','Color',cmap(thisColor,:));
    hold on
    plot(theseTimes,thisTraj,'.','Color',cmap(thisColor,:));
    
    
%     plot(theseTimes,thisTraj,'.-','Color',cmap(thisColor,:));
%     semilogy(theseTimes,thisTraj,'.-','Color',cmap(thisColor,:));
    hold on
    
    
  end
 
  
  %deterministic growth
    detTime = -100:10:300;
    threshRat = thresh/(1-thresh);
    detTraj = threshRat*exp(detTime*s(i));
    detTraj  = detTraj./(1+detTraj );
    plot(detTime,detTraj,'--','LineWidth',3,'Color',[.75 .75 .75])  
    
    
%     detTraj = exp(log(thresh/(1-thresh)) + detTime*s(i));
%     plot(detTime,detTraj./(detTraj+1),'c.','MarkerSize',20)
    
  
  
%   theseTimes = 0:50:750;
%   thisInoc = mean(inocFreqBestFit(pickThese));
%   determGrowth = (thisInoc /(1-thisInoc ))*exp(s(i)*theseTimes );
%   determGrowth = determGrowth./(1+determGrowth);
%   plot(theseTimes , determGrowth,'k--','LineWidth',3)
  
  %cosmetics
  axis([thisXaxis 0 1]);
  set(gca,'FontSize',14)
  caxis([min(colorScale) max(colorScale)])
  
%    colorbar
end

  
  for i=1:4
      subplot(2,2,i)
      ylim([0.01 1])
      
      if i==3
          xlim([-100 650])
      end
  end

clear this* these* i j pick* s
  
%%
for i=1:4
    subplot(3,4,i)
    colorbar('off')
end

%%

titText= {'2.8','3.9','5.0','7.3'};
for i=1:4
    subplot(2,2,i)
    if i==3 || i==4
        xlabel('generations')
    else xlabel('');
    end
    if i==1 || i==3
        ylabel('frequency of seeded lineage')
    else ylabel('')
    end
        
    % title(['s_{seed} = 2.8%, exponential DFE mean ' num2str(DFE_params(2)) ' mutation rate 10^{' num2str(DFE_params(1)) '}'])  
%     title('experimental data trajectories for s_{seed} = 2.8%')
    title(['fitness = ' titText{i} '%']) 
    colorbar('off')
end



%% plot all trajs

s = fitClasses;
thresh = .1;

%cosmetics
cmap = colormap; 
cmap = get(gcf,'Colormap');
close(gcf);
colorScale = linspace(...
        -3,...
        -2,...
        size(cmap,1));

%for each fitness class
figure
for i=3
  pickThese = find(trajFitClass  == fitClasses(i));
  pickThese = pickThese(randperm(length(pickThese)));

  %time range
  thisXaxis = [0 0];
   
  
  for j=1:length(pickThese)
      
    theseTimes =  trajs{pickThese(j)}(1,:)*10; 
    thisInocYR =  trajs{pickThese(j)}(2,1);    
    thisInocFreq = inocFreqBestFit(pickThese(j));
    
    
     thisTraj = trajs{pickThese(j)}(2,:)./(trajs{pickThese(j)}(2,:)+1); 

    
    timeInd2 = find(thisTraj>thresh,1,'first');
    timeInd1 = find(thisTraj(1:timeInd2)<thresh ,1,'last');
    timeInds  =  [timeInd1 timeInd2 ]; clear timeInd1 timeInd2
    theseY = thisTraj(timeInds);
    theseX = theseTimes(timeInds);
    if length(timeInds)==2
        deltaT = ( theseX(1)*(theseY(2)-thresh) + theseX(2)*(thresh - theseY(1)) )...
            /( theseY(2) - theseY(1));
    else
        thisT = (1/s(i))*log((thresh/(1-thresh))/inocFreqBestFit(pickThese(j)));
    end

    theseTimes = theseTimes-deltaT;    
    thisXaxis = [min(thisXaxis(1),theseTimes(1)) max(thisXaxis(2),theseTimes(end))];
    
    
    thisColor = logical(hist(log10(thisInocYR),colorScale));
       
    plot(theseTimes,thisTraj,'.-','Color',cmap(thisColor,:));
    hold on
    
  end
  
    
  %deterministic growth
    detTime = -100:10:300;
    detTraj = exp(log(thresh/(1-thresh)) + detTime*s(i));
%     plot(detTime,detTraj./(detTraj+1),'c.','MarkerSize',20)
    plot(detTime,detTraj./(detTraj+1),'k--','LineWidth',3)  

  
  %cosmetics
  axis([thisXaxis 0 1]);
  set(gca,'FontSize',14)
  caxis([min(colorScale) max(colorScale)])
  
%    colorbar
end
xlim([-200 650])
clear this* these* i j pick* s








%%
% %simulate pfix for three (s,inocs)
% 
% % theseFits = [2.8 5.2 7.3]*10^-2;
% % theseInocs = [.6 .55 1.6]*10^-2;
% %cd 'C:\Users\Genya\Documents\Eclipse workspace\BMA sims Matlab\likeParamEstimate\likelyParamEstimate_deltaDFE'
% 
% trials = 10;
% 
% %trajs = cell{size(theseFits)};
% 
% %pick = 251;
% 
% Nb= 2*10^4;
% gPX = 10;
% U=10^-4;
% alpha = .03;
% s = .04*(2*rand(1,trials)+.5);
% measTimes = 0:60:1000;
% inoc = .01*(rand(1,trials) + .25);
% 
% %exponential MLEs
% % U=10^-5.5;
% % alpha = .026;
% %exponential MLE for s=3% fitness class
% % U = 10^-4.6;
% % alpha = .011;
% 
% %dirac delta MLEs
% % U=10^-3;
% % alpha = .04;
% 
% 
% % s = .028*ones(1,trials);
% % inoc = .006*ones(1,trials);
% 
% %select some random meas times & inocs
% % pick = find(trajFitns < .045 & trajFitns > .035);
% % pick = find(trajFitns == .028);
% % pickRand = randperm(length(pick));
% % pickRand = pick(pickRand);
% % measTimes = cell(trials,1);
% % for i=1:trials
% % %     inoc(i) = trajs{pickRand(i)}(2,1)./(trajs{pickRand(i)}(2,1)+1);
% %     measTimes{i} = trajs{pickRand(i)}(1,:)*10;
% % end
% 
% %s = trajFitns(258);
% %measTimes = trajs{258}(1,:)*10;
% %inoc = trajs{258}(2,1)./(trajs{258}(2,1)+1);
% 
% % simTrajs = zeros(trials,length(measTimes));
% simTrajs =cell(trials,1);
% % simFitTrajs = zeros(trials,length(measTimes));
% % timeTrajs = zeros(1,trials);
% 
% for i=1:trials
% 
%     %[simTrajs(i,:) simFitTrajs(i,:)] = simTraj7_expDFE(s,inoc,measTimes,Nb,gPX,alpha,U);
%     %simTrajs{i} = simTraj8_expDFE(s(i),inoc(i),measTimes{i},Nb,gPX,alpha,U);
% %     timeTrajs(i) = toc;
%      simTrajs{i} =  simTraj8(s(i),inoc(i),measTimes,Nb,gPX,U,alpha,'dir');
%    % simTrajs{i} =  simTraj7_expDFE(s(i),inoc(i),measTimes,Nb,gPX,alpha,U);
% 
% end
% %simTrajs
% % mean(timeTrajs)
% %%
% 
% colorList = [
% 	0.00  0.00  1.00
% 	0.00  0.50  0.00
% 	1.00  0.00  0.00
% 	0.00  0.75  0.75
% 	0.75  0.00  0.75
% %	0.75  0.75  0.00
% 	0.1  0.1  0.1
% 	0.75  0.25  0.25
% 	0.95  0.95  0.00
% %	0.25  0.25  0.75
% %	0.75  0.75  0.75
% 	0.00  1.00  0.00
% 	0.76  0.57  0.17
% 	0.54  0.63  0.22
% 	0.34  0.57  0.92
% 	1.00  0.10  0.60
% %	0.88  0.75  0.73
% 	0.10  0.49  0.47
% %	0.66  0.34  0.65
% 	0.99  0.41  0.23
% ]; 
% %%
% %all meastimes
% allMeasTimes = [];
% for i=1:trials
%     allMeasTimes = unique([allMeasTimes measTimes]);
% end
% %deterministic growth
% inocRat = (1-mean(inoc))/mean(inoc);
% expGrow = 1./(1+ inocRat*exp(-mean(s)*allMeasTimes));
% 
% figure
% colors = varycolor(size(simTrajs,1));
% for i = 1:size(simTrajs,1)
%     plot(measTimes,simTrajs{i},'.-','Color',colorList(mod(i,size(colorList,1))+1,:),...
%         'LineWidth',2,'MarkerSize',15)
%     hold all
% end
% plot(allMeasTimes,expGrow,'k--','LineWidth',2.5)
% hold off
% xlim([0 800])
% set(gca,'FontSize',14)
% %%
% 
% %plot(simFitTrajs')
% 
% %get pfix
% % fixers = zeros(trials,1);
% % for i=1:trials
% %     if sum(trajs(i,:) == 1)>0
% %         fixers(i) = 1;
% %     end
% % end
% % sum(fixers);
% 
% % for i=1:size(trajs,1)
% %    if trajs(i,8) > .5
% %        trajs(i,9) = 1;
% %    else
% %        trajs(i,9) = 0;
% %    end
% % end
% %%
% title('')
% xlabel('')
% ylabel('')
% 
% %%
% figure
% plot(measTimes,simTrajs,'.-','MarkerSize',15)
% % hold on
% % plot(measTimes,expGrow,'k--','LineWidth',3)
% % hold off
% xlim([0 1000])
% fs = 13;
% xlabel('generations','FontSize',fs)
% ylabel('frequency of seeded lineage','FontSize',fs)
% set(gca,'FontSize',fs)
% title(['\rho(s) = e^{-s/' num2str(alpha) '}, U_b=10^{' num2str(log10(U)) '}' ],'FontSize',fs)
% %title(['\rho(s) = \delta(' num2str(alpha) '-s), U_b=' num2str(U) ],'FontSize',fs)
%     
%scripts to test returnGenMutFunct() for different DFEs
clear
addpath('simulation code')
%% exponential DFE

DFE_form  = 'exp';
DFE_params = [-1 .05];

[genMuts delS] = returnGenMutFunct(DFE_form,DFE_params);

realPDF = @(x,mu) (1/mu)*exp(-x/mu);

n = 10^5;

[mutsPerPop fitSteps] = genMuts(n);

[thisPDF bins] = hist(fitSteps*delS,round(mutsPerPop/100));

dBin = bins(2)-bins(1);

thisPDF = (thisPDF/mutsPerPop)/dBin;

stairs(bins,thisPDF)
hold on
plot(bins,realPDF(bins,DFE_params(2)) ,'r--')
xlabel('finess')
ylabel('probability')
title(['U \rho(s) = U\cdot(1/\mu)\cdot e^{-s/\mu}, \mu=' num2str(DFE_params(2)) ', U = 10^{' num2str(DFE_params(1)) '}'])
legend('from genMutFunct()','analytical')

%% dirac delta

DFE_form  = 'dir';
DFE_params = [-1 .05];

[genMuts delS] = returnGenMutFunct(DFE_form,DFE_params);

n = 100;

[mutsPerPop fitSteps] = genMuts(n);

mutsPerPop/(n*10^DFE_params(1))

std(fitSteps)
[mean(fitSteps) DFE_params(2)/delS ]


%% truncated exponential DFE
 
DFE_form  = 'truncExp';
DFE_params = [-2 .05 .02 .1];

[genMuts delS] = returnGenMutFunct(DFE_form,DFE_params);

realPDF = @(x,mu,xMin,xMax) (heaviside(x-xMin)/mu).*exp(-x/mu)...
    /(exp(-xMin/mu)-exp(-xMax/mu));

n = 10^6;

[mutsPerPop fitSteps] = genMuts(n);

[thisPDF bins] = hist(fitSteps*delS,round(mutsPerPop/100));

dBin = bins(2)-bins(1);

thisPDF = (thisPDF/mutsPerPop)/dBin;

stairs(bins,thisPDF)
hold on
plot(bins, realPDF(bins,DFE_params(2),DFE_params(3),DFE_params(4)),'r--')
xlim([0 DFE_params(4)*1.2])

xlabel('finess')
ylabel('probability')
title(['U \rho(s) = U \cdot (1/\mu)\cdot e^{-s/\mu} truncated at ' num2str(DFE_params(3)) ' from the left and ' num2str(DFE_params(4)) ' from the right, \mu=' num2str(DFE_params(2)) ', U = 10^{' num2str(DFE_params(1)) '}'])
legend('from genMutFunct()','analytical')

%% uniform DFE

DFE_form = 'unif';
DFE_params = [-1 .03];
[genMuts delS] = returnGenMutFunct(DFE_form,DFE_params);

realPDF = @(x,mu) heaviside(2*mu-x)/(2*mu);

n=10^7;

[mutsPerPop fitSteps] = genMuts(n);

[thisPDF bins] = hist(fitSteps*delS,round(mutsPerPop/10^4));

dBin =  bins(2)-bins(1);

thisPDF = (thisPDF/mutsPerPop)/dBin;

stairs(bins,thisPDF)
hold on
plot(bins,realPDF(bins,DFE_params(2)) ,'r--')
hold off
xlim([0 DFE_params(2)*2.2])
xlabel('fitness')
ylabel('probability')
title(['U \rho(s) = U \cdot \Theta_{Heaviside}( 2\cdot\mu - s ), \mu = ' num2str(DFE_params(2)) ', U = 10^{' num2str(DFE_params(1)) '}'])
legend('from genMutFunct()','analytical')


%% piecewise DFE

DFE_form  = 'piecewise';
DFE_params = [-1 -1.25 -1.5 -1.75 -2]; %MUTATION RATE FOR EACH INTERVAL
[genMuts delS] = returnGenMutFunct(DFE_form,DFE_params);

    %mutation rates for mutations in intervals defined by 
    DFE_parts = [.028; .038; .05; .073; .15];
    
    dS = .002; %bin width for histogram
    
    DFEbinEdgs = [0:dS:DFE_parts(end)];
    DFEbinCents = [dS/2:dS:DFEbinEdgs(end)-dS/2];
    sRange = DFEbinEdgs(end) - DFEbinEdgs(1);
    DFE_partWidths = DFE_parts - [0; DFE_parts(1:end-1)];
    totMutRate = sum(10.^DFE_params);
    DFE_weights = 10.^DFE_params/totMutRate;  
    
n = 10^5;
[mutsPerPop fitSteps] = genMuts(n);
[mutsBinned bins] = histc(fitSteps*delS,DFEbinEdgs);
mutsBinned = mutsBinned(1:end-1);
probDens = mutsBinned/(n*dS);

predHist = zeros(length(DFEbinCents),1);
for i=1:length(DFEbinCents)
    pick = find(DFEbinCents(i) < DFE_parts,1);
    predHist(i) = n*totMutRate*DFE_weights(pick)*(DFEbinEdgs(i+1)-DFEbinEdgs(i))/DFE_partWidths(pick);
end

bar(DFEbinCents,mutsBinned)
hold on
stairs(DFEbinCents,predHist,'r','LineWidth',3)
% hold all
% for i=1:length(DFE_parts)
%    vline(DFE_parts(i),'g--')
% end

xlim([DFEbinEdgs(1) DFEbinEdgs(end-1)*2])
xlabel('fitness')
ylabel('number of mutations')
legend('returnGenMutFunct()','analytical')

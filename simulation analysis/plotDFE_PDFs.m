%plot shape inferred most-likely DFEs

figure
x = .02:.0001:.073;

%cosmetics
lineWid = 3;
fs = 35;

%exponential
expMu = 0.0085;
expU = 10^-4;
expPDF = @(x,mu) (1/mu)*exp(-x/mu);
expY = expU*expPDF(x,expMu);
 plot(x,expY,'r-',...
     'LineWidth',lineWid )
hold all

%uniform
unifMu = 0.032;
unifU = 10^-5.69;
unifPDF = @(x,mu) heaviside(2*mu-x)/(2*mu);
 plot(x,unifU*unifPDF(x,unifMu),'b-',...
     'LineWidth',lineWid )

%delta
dirMu = 0.0585;
dirU = 10^-6.29;
bin = .0001;
height = 1/bin;
dirX = [dirMu - bin/2, dirMu - bin/2, dirMu + bin/2, dirMu + bin/2];
dirY = [0 height height 0];
plot(dirX,dirY*dirU,'c-',...
     'LineWidth',lineWid )

%more cosmetics
xlim([0 .08])
ylim([0 max(expY)*1.1])
ylabel('U\rho(s)')
xlabel('mean of DFE, \rho(s)')

formatFig(gcf,gca,fs)
legend('exponential','uniform','\delta-function')






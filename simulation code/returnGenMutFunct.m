function [genMuts delS] = returnGenMutFunct(DFE_form,DFE_params)
%function [genMuts delS] = returnGenMutFunct(DFE_form,DFE_params)
%
% genMuts(n) = [mutsPerPop mutSteps] where
% mutsPerPop = poissrand(sum(U)*n) and 
% mutSteps = an array of length sum(mutsPerPop) containing random numbers 
% drawn from a distribution specified by DFE_form and rounded to 
% 1 + the nearest multiple of 1/delS. 
%
% delS = 10^-4, except for dirac delta DFE.
%
% DFE_params has different meanings depending on DFE_form
%
% if DFE_form = 'truncExp', 
%    then DFE_params = [mutation rate; mean of DFE; cut-off]
%
% if DFE_form = 'exp' or 'dir' 
%    then DFE_params = [mutation rate; mean of DFE]
%
% if DFE_form = 'piecewise', 
%   the DFE is divided into uniform density intervals:
%       [0 to DFE_parts(1) ), [DFE_parts(2) to DFE_parts(3)), etc. 
%       where DFE_parts = [.028, .038, .05, .073]
%   and DFE_params = [mutation rates in each of these intervals],
%
% if DFE_form = 'gamma', DFE_params = [mutation rate; shape param; scale param]

    %parse inputs (defines/computes vars used by DFE functs)
    if strcmp(DFE_form,'exp')
        delS = 10^-4;
        U = 10^DFE_params(1);
        alpha = DFE_params(2);
        if length(DFE_params) ~= 2
            error('for exponential DFE, not right number of params specified');
        end
    elseif strcmp(DFE_form,'dir')
        delS = 10^-4;
        U = 10^DFE_params(1);
        delStep = round(DFE_params(2)/delS);
        if length(DFE_params) ~= 2
            error('for dirac delta DFE, not right number of params specified');
        end
    elseif strcmp(DFE_form,'truncExp')
        delS = 10^-4;
        U = 10^DFE_params(1);
        alpha = DFE_params(2);
        truncL = DFE_params(3);
        truncR = DFE_params(4);
        if length(DFE_params) ~= 4
            error('for truncated exponential DFE, not right number of params specified');
        end
    elseif strcmp(DFE_form,'piecewise')
        DFE_parts = [.028; .038; .05; .073; .15];
        if length(DFE_params) ~= length(DFE_parts)
            error('for piecewise DFE, too many or few mutation rates specified');
        end
        
        U = 10.^DFE_params(:);
        delS = 10^-4;
        totU = sum(U);
        
        fitClasses = [0; DFE_parts]/delS;
        delClass = (fitClasses(2:end) - fitClasses(1:end-1));

        DFE_weights = U/totU;

        dfeCDF = [0; cumsum(DFE_weights)];
        delDFEcdf = (dfeCDF(2:end) - dfeCDF(1:end-1));
    elseif strcmp(DFE_form,'gamma')
        if length(DFE_params) ~= 3
            error('for gamma DFE, not right number of params specified');
        end        
        
        U = 10^DFE_params(1); %mutation rate
        gamK = DFE_params(2); %shape parameter
        gamT = DFE_params(3); %scale parameter
    
    elseif strcmp(DFE_form,'unif')
               
        U = 10^DFE_params(1); %mutation rate
        alpha = DFE_params(2); %mean
        
        delS = 10^-4;
    else
        error('input string ''DFE_form'' is unrecognized');
    end    
       
    %define DFE functions
    function [mutsPerPop fitSteps] = expDFE(n)
        mutsPerPop = poissrnd(U*n);
        fitSteps = round(exprnd(alpha,[sum(mutsPerPop) 1])/delS);
    end

    function [mutsPerPop fitSteps] = dirDFE(n)
        mutsPerPop = poissrnd(U*n);
        fitSteps = ones(sum(mutsPerPop),1)*delStep;
    end

    function [mutsPerPop fitSteps] = truncExpDFE(n)
        %an exponential with the same mean and mutation rate but truncted
        
        %weight the mutation rate by the probability included within the
        %trunction bounds
        mutsPerPop = poissrnd( U*(exp(-truncL/alpha) - exp(-truncR/alpha))*n);               
        
        
        invTruncExpCDF = @(F,mu,tMin,tMax) -mu*log(exp(-tMin/mu)-F*(exp(-tMin/mu)-exp(-tMax/mu)));
        
        fitSteps = round(...
            invTruncExpCDF(rand(sum(mutsPerPop),1),alpha,truncL,truncR) ...
            /delS);
    end

    function [mutsPerPop fitSteps] = pieceDFE(n)
        
       mutsPerPop = poissrnd(totU*n);
       totMuts = sum(mutsPerPop);
        
       fitSteps = rand(totMuts,1);
       k = zeros(totMuts,1);
        
        i = 2;
        while find(k==0,1)
              k( (fitSteps >= dfeCDF(i-1)) & (fitSteps < dfeCDF(i)) )  = i-1;
              i = i + 1;
        end
        
        fitSteps = floor(  ...
            fitClasses(k) + delClass(k).*(fitSteps - dfeCDF(k))./delDFEcdf(k)...
                );
    end

    function [mutsPerPop fitSteps] = gamDFE(n)
        mutsPerPop = poissrnd(U*n);
        fitSteps = round(gamrnd(gamK,gamT,[sum(mutsPerPop) 1])/delS);
    end

    function [mutsPerPop fitSteps] = unifDFE(n)
        mutsPerPop = poissrnd(U*n);
        fitSteps = ceil(rand([sum(mutsPerPop) 1])*2*alpha/delS);
    end

    %return appropriate DFE function
    if strcmp(DFE_form,'exp')
        genMuts = @expDFE;
    elseif strcmp(DFE_form,'dir')
        genMuts = @dirDFE;
    elseif strcmp(DFE_form,'truncExp')
        genMuts = @truncExpDFE;
    elseif strcmp(DFE_form,'piecewise')
        genMuts = @pieceDFE;
    elseif strcmp(DFE_form,'gamma')
        genMuts = @gamDFE;        
    elseif strcmp(DFE_form,'unif')
        genMuts = @unifDFE;                
    else
        error('input string ''DFE_form'' is undefined');
    end

end
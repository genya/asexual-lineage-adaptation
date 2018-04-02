function [xtremaBin peakFreqSdownBin] = classifyTraj(fullTraj,allMeasTimes,xtremaBins,peakFreqSdownBins,subset)
%[xtremaBin peakFreqSdownBin] = classifyTraj(fullTraj, ...
%                                            allMeasTimes, ...
%                                            xtremaBins, ...
%                                            peakFreqSdownBins, ...
%                                            subset);
%
% fullTraj is assumed to be frequencies not ratios
%
% [xtremaBin peakFreqSdownBin] are scalars
%
% peakFreqSdownBin is assumed to be of the form
% [first peak freq bin centers s_down bin centers]
%
% Code assumes that no consecutive measurement times are more than 110
% generations apart. May also break if the first peak in the trajectory 
% consists of two consecutive measurements that are exactly the same. 


%reduce the trajectory to a desired subset of time points
if nargin < 5
    measTimes = allMeasTimes;
    traj = fullTraj;
else
    traj = fullTraj(logical(subset));
    measTimes = allMeasTimes(logical(subset));
end
if nargin < 4
    peakFreqSdownBins = [];
end
if nargin < 3
    error('not enough input arguments')
end
     
%determine # of extema
fluctThresh = .04; %hard-coded param
[maxtab mintab] = peakdet(traj, fluctThresh);
numXtrema = size(maxtab,1) + size(mintab,1);
[~, xtremaBin] = min(abs(xtremaBins - numXtrema));

%determine (peak freq, s_down) bin
if isempty(peakFreqSdownBins)
    peakFreqSdownBin = [];
elseif numXtrema == 0 
   %if no peak, sort trajs into those with high final frequency 
   %(essentially monotonically increasing) and those with low final frequ 
   %("no rise" trajectories)
    
   %which bins have s_down undefined
   pick = find(isnan(peakFreqSdownBins(:,2)));
   
   %among these, which has closest final freq
   [~, peakFreqSdownBin] =  min( abs( ...
       peakFreqSdownBins(pick,1)-fullTraj(end)   ));
   
   %find the index of that bin
   peakFreqSdownBin = pick(peakFreqSdownBin);
   
   if isempty(peakFreqSdownBin)
       peakFreqSdownBin = NaN;
   end
else

    %select peak & downslope time pts
    peakTime = measTimes(maxtab(1,1));
        %use the measurement ~60 generations after the peak to determine s_down
        %this is an AD HOC fix to address the statistical bias 
        %introduced by the fact that I stopped measuring
        %trajectories after they resolved. Simulated trajectories therefore
        %sometimes have added measurement times, corresponding to when they would 
        %have been measured had the seeded lineage not yet fixed or gone extinct.
    [~, pickNext] = min( abs( (allMeasTimes-peakTime)-60) );
    postPeakTime = allMeasTimes(pickNext);
    slopeMeasTimes = [peakTime postPeakTime];
    slopeTraj = [ traj(maxtab(1,1)) fullTraj(pickNext)];                 
        
    %compute s_down
    slopeTraj = slopeTraj./(1-slopeTraj);
    sdown = abs(... %sdown is positive by convention
                        log(slopeTraj(2)/slopeTraj(1))...
                    /( slopeMeasTimes(2) - slopeMeasTimes(1)));    
                
   %round s_down to nearest bin
   [~, thisSdownBin] = min(abs(peakFreqSdownBins(:,2)-sdown));
    thisSdownBin = peakFreqSdownBins(thisSdownBin,2);
   
   %identify closest peak bin 
   thesePeakBins = find(thisSdownBin == peakFreqSdownBins(:,2));
   [~, thisPeakBin] = min(abs(...
       peakFreqSdownBins(thesePeakBins,1)-traj(maxtab(1,1))   ));
    thisPeakBin = peakFreqSdownBins(thesePeakBins(thisPeakBin),1);
        
    %find bin in joint dist of peak height & s_down
    peakFreqSdownBin = find( ...
        peakFreqSdownBins(:,1) == thisPeakBin...
        & peakFreqSdownBins(:,2) == thisSdownBin ,1);
   
end


    function [maxtab, mintab]=peakdet(v, delta)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    % [Slightly modified by EMF, 9/24/2012]

    maxtab = [];
    mintab = [];

    v = v(:); % Just in case this wasn't a proper vector
    x = (1:length(v))';

    if (length(delta(:)))>1
      error('Input argument DELTA must be a scalar');
    end

    if delta <= 0
      error('Input argument DELTA must be positive');
    end

    mn = Inf; mx = -Inf;
    mnpos = NaN; mxpos = NaN;

    lookformax = 1;

    for i=1:length(v)
      this = v(i);
      if this > mx, mx = this; mxpos = x(i); end
      if this < mn, mn = this; mnpos = x(i); end

      if lookformax
        if this < mx-delta
          maxtab = [maxtab ; mxpos mx];
          mn = this; mnpos = x(i);
          lookformax = 0;
        end  
      else
        if this > mn+delta
          mintab = [mintab ; mnpos mn];
          mx = this; mxpos = x(i);
          lookformax = 1;
        end
      end
    end
    end
    
 end
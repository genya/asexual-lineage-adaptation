function [traj] = adapRate(measTimes,genMuts,delS)
%function [traj] = adapRate(measTimes,genMuts,delS)
%
% returns mean fitness at measTimes of an initially clonal, s=0 pop.

%hard-coded params
Nb = 10^4;
gPX = 10;
smallThresh = 10; 
    %dilution of fitness classes with more cells than this use gaussian
    %approx to poisson distribution

%initialize output
traj = zeros(length(measTimes),1)-1;
    %traj(i) = -1 => traj(i) not yet recorded

Nold = 1;
Nmax = 3*10^8;
sizes = sparse(1,... %each row is a fitness class, fitness = (row index - 1) * delS
        1,... %each column is a lineage
        Nold); %starting pop sizes
t= 0; 
while Nold <= Nmax

    %calulate growth (deterministic)
    [fitInds, ~, cnts] = find(sizes);
    Nold = sum(cnts);

    tStep = min(2,max(log2(Nmax/Nold),0));  %number of doublings in this time step

    Nnew = (2^tStep)*cnts;%everyone grows neutral %.*exp(fits*tStep);

    %draw mutants    
    [muts fitSteps] = genMuts(Nnew);
    totMuts = sum(muts);    

    if totMuts > 0

        %distribute mutations
            %for each fitStep, need to add the index of its source fitness
            %class and assign it to its source lineage   
        [mutInds, ~, muts] = find(muts);

        x = 1;
        for i=1:length(mutInds)
            pick = x:x+muts(i)-1;
            fitSteps(pick) = fitSteps(pick) + fitInds(mutInds(i)); 
            x = x + muts(i);
        end    

        %combine old guys, not-mutated offspring, and mutants 
        sizes = sparse(...
            [fitInds; fitSteps],...
            ones(totMuts+length(fitInds),1),...
            [Nnew; ones(totMuts,1)]); 
    else
        %update pop
        sizes = sparse(fitInds,ones(length(fitInds),1),Nnew); 
    end

    %increment time
    t = t + tStep;   

end    
    
%subsample
[fitInds, ~, cnts] = find(sizes);
cnts = poissrnd(cnts*Nb/sum(cnts));

sizes = sparse(fitInds,...
            ones(length(cnts),1),...
            cnts);     

%t = 0
if ~isempty(find(measTimes==0,1))
    traj(measTimes==0) = 0;
    measTimes(measTimes==0) = -1;
end
nextTimePt= find(measTimes>0,1);

%run simulation
t=0;
while t<= measTimes(end)

    %deterministic growth & dilution
    [fitInds, linInds, cnts] = find(sizes);    
    
    %record trajectory
    if t >= measTimes(nextTimePt)    
        traj(nextTimePt) = sum(cnts.*(fitInds-1)*delS)/sum(cnts);
        measTimes(nextTimePt) = -1; %mark current time pt as measured
        nextTimePt = find(measTimes>0,1); %get next time pt
    end         
    
   
%     Z = cnts.*exp(gPX*(fitInds-1)*delS); 
%     Nnew = Nb*Z/sum(Z); %expected fitness class sizes after dilution
    Nnew = Nb*(cnts.*exp(gPX*(fitInds-1)*delS))...
          /sum(cnts.*exp(gPX*(fitInds-1)*delS));

    %stochastic dilution
    pickSmall = Nnew < smallThresh;
    
    Nnew( pickSmall) = poissrnd(Nnew(pickSmall));
    Nnew( ~pickSmall) = randn(size(Nnew(~pickSmall))).*sqrt(Nnew(~pickSmall)) + Nnew(~pickSmall);
    Nnew(Nnew < 0) = 0;

    %draw mutants    
    [muts fitSteps] = genMuts(Nnew*gPX);
    totMuts = sum(muts);
    
    %draw mutational effects
    if  totMuts > 0

        Nnew = Nnew-muts; Nnew(Nnew < 0) = 0;

        %distribute mutations
            %for each fitStep, need to add the index of its source fitness
            %class and assign it to its source lineage   
        mutLin = zeros(totMuts,1);
        [mutInds, ~, muts] = find(muts);

        x = 1;
        for i=1:length(mutInds)
            pick = x:x+muts(i)-1;
            fitSteps(pick) = fitSteps(pick) + fitInds(mutInds(i)); 
            mutLin(pick) = linInds(mutInds(i));
            x = x + muts(i);
        end    
        
        %update pop with mutants 
        sizes = sparse(...
            [fitInds; fitSteps],...
            [linInds; mutLin],...
            [Nnew; ones(totMuts,1)]); 
    else
        %update pop
        sizes = sparse(fitInds,linInds,Nnew); 
    end
    
    %increment time
     t = t + gPX;             
end
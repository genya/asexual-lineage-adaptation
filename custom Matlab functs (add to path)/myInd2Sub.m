function subscripts = myInd2Sub(matSize,linInd)
%subscripts = myInd2Sub(matSize,linInd)
%
% returns same values as ind2sub but in a different form. Does not accept
% non-scalar linInd.
%
% If  ind2sub(SIZ,IND)  = [I1,I2,I3,...,In] 
% then myInd2Sub(SIZ,IND)  = I
%  where I = [I1,I2,I3,...,In] and I1,I2,I3,...,In are scalars.
% 

matSize = matSize(:);

if length(matSize) == 1
    subscripts = linInd;
else

    subscripts = zeros(length(matSize),1);

    dimProd = cumprod(matSize);
    dimProd = [0; dimProd(1:end-1)];

    thisInd = linInd;
    for i=0:length(subscripts)-1

        if dimProd(end-i) > 0
            
            if thisInd <= dimProd(end-i)
                subscripts(end-i) = 1;
            else
                if mod(thisInd,dimProd(end-i)) == 0
                    subscripts(end-i) = thisInd/dimProd(end-i);        
                else               
                    subscripts(end-i) = floor(thisInd/dimProd(end-i)) + 1;    
                end
                thisInd = thisInd - (subscripts(end-i)-1)*dimProd(end-i);                                        
            end
        else
            
            subscripts(end-i) = max(thisInd,1);
        end
    end
   
end

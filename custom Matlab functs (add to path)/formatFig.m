function formatFig(figHan,axHan,varargin)
%formatFig(figHan,axHan,fontSize)
% 
% fontSize is optional input. default value is 13

if isempty(varargin)
    fs = 13;
else
    fs = varargin{1};
end

set(axHan,'FontSize',fs)
set(get(axHan,'YLabel'),'FontSize',fs)
set(get(axHan,'XLabel'),'FontSize',fs)
set(get(axHan,'title'),'FontSize',fs)

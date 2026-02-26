

SPN = 'F:\S1-1Kx1Kcore_one_second\'

%% Parse file names

oldNames = dir([SPN '*.tif']);
clear oldSec
for i = 1:length(oldNames)
   nam = oldNames(i).name; 
   oldSec(i,:) = sscanf(nam,'z%d_W%s',1);
end

%%
%secondOrder = [];
for i = 1:size(secondOrder,1)
    targ = find(oldSec == secondOrder(i,1));
    if isempty(targ)
        'no file found'
    elseif length(targ)>1
        'more than one targ found'
        break
    else
        oldNam = oldNames(targ).name;
        newNam = sprintf('%sz%05.0f%s',SPN,secondOrder(i,4),oldNam(7:end))
       movefile([SPN oldNam ],newNam)
    end
end
    

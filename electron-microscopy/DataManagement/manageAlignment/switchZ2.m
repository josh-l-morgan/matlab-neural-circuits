

%SPN = 'D:\LGNs1\PSC_alignments\Nov_4_1Kx1Kcore - Copy\'

%% Parse file names

oldNames = dir([SPN '*.tif']);
clear oldSec
tags = {'_','_','_'}

for i = 1:length(oldNames)
   nam = oldNames(i).name; 
   %oldSec(i,:) = sscanf(nam,'%*s_%d_%d_%d%s',1);
    oldSec(i,:) = parseNum(nam,tags);
end

%%
%secondOrder = [];
for i = 1:size(secondOrder,1)
   %targ = find(oldSec == secondOrder(i,1));
    targ = find((oldSec(:,2) == secondOrder(i,2)) & (oldSec(:,3) == secondOrder(i,3)))
    if isempty(targ)
        'no file found'
    elseif length(targ)>1
        'more than one targ found'
        break
    else
        oldNam = oldNames(targ).name
        newNam = sprintf('%s_%05.0f_%02.0f_%03.0f_core.tif',oldNam(1:6),secondOrder(i,4),secondOrder(i,2),secondOrder(i,3))
      try
        movefile([SPN oldNam ],[SPN newNam])
      catch err
          err
      end
    end
end
    

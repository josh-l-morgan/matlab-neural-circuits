
%
%clear all
SPN = GetMyDir;
dSPN = dir([SPN '*.tif'])
TPN =[SPN(1:end-1) '_reorder\'];
mkdir(TPN)

%%
dir([SPN '*.tif']);
clear w s
for i = 1:length(dSPN)
    nam = dSPN(i).name;
    und = regexp(nam,'_');
    w(i) = str2num(nam(und(1)+1:und(2)-1));
    s(i) = str2num(nam(und(2)+1:end-4));
end
%%
o = 1:length(s);

for i = 1:size(reorder,1);
    neworder = reorder(i,:);
    neworder = neworder(neworder>0);
    if ~isempty(neworder)
       
    waf = neworder(1);
    sec = neworder(2:end);
    clear secO
    for si = 1:length(sec)
        foundO =  find(((w==waf) + (s==sec(si)))==2);
        if isempty(foundO)
            secO(si) = 0;
        else
            secO(si) = foundO;
        end
    end
    secO = secO(secO>0);
    if length(secO)>1
        sortSecO = sort(secO);
       o(secO) = sortSecO; 
    end
    end
    
end

allOrders = [w' s' o']
%%
%{
imageOrder = 1:length(o);
imageOrder = imageOrder(o); %index is image order, value is position within dSPN

for i = 1:length(imageOrder)
  oldName =  dSPN(imageOrder(i)).name; 
  newName = sprintf('%05.0f_%03.0f_%03.0f.tif',i,w(imageOrder(i)),s(imageOrder(i)));
  if ~exist([TPN newName],'file')
      try
        copyfile([SPN oldName],[TPN newName]);
      catch err
          err
      end
  i
  end
end

%}

%%
for i = 1:max(o)
    whosNext = find(o == i);
    if ~isempty(whosNext)
        if length(whosNext)>1
            disp('too many')
        else
            oldName =  dSPN(whosNext).name;
            
            newName = sprintf('%05.0f_%03.0f_%03.0f.tif',i,w(whosNext),s(whosNext));
            
            if ~exist([TPN newName],'file')
                try
                    copyfile([SPN oldName],[TPN newName]);
                catch err
                    err
                end %try copy
                i
            end % if no target file
        end %if one target
        
    end % if some target
end %run all 
    


%%
% 
% 
%     
%     tifInd = str2num(nam(end-7:end-4));
%     newname = sprintf('%s%04.0f.tif',basename,tifInd);
%     %     copyfile([SPN nam],[TPN newname]);
%     %     if (i <L) && ( i>1)
%     %     newInd = L+ L-i;
%     %     newername = sprintf('%s%04.0f.tif',basename,newInd);
%     %     copyfile([SPN nam],[TPN newername]);
%     
% end
% end
% 

return
SPN = 'C:\myData\LGNs1\LGNS1_alignmentData\s8uploadSamp\S8uploadingSubSampMiddle256\';
TPN = 'C:\myData\LGNs1\LGNS1_alignmentData\s8uploadSamp\S8uploadingSubSampMiddle256_reordered\'
if ~exist(TPN,'dir'),mkdir(TPN),end

load('.\newOrder.mat')

%%
tifs = dir([SPN '*.tif']);

for i = 1:length(tifs)
    
    nam = tifs(i).name;
    W = regexp(nam,'W');
    S = regexp(nam,'_s');
    dotTif = regexp(nam,'.tif');
    w = str2num(nam(W(end) + 1:S(end) - 1));
    s = str2num(nam(S(end) + 2:dotTif(end)-1));
    secID(i) = str2num(sprintf('%02.0f%04.0f',w,s));
    
    
end

%%
secOrder = secID;
for i = 1:length(newOrder)
    w = newOrder(i,1);
    secs = newOrder(i,2:end);
    secs = secs(secs>0);
%     if (w == 32) & length(secs)>4
%         pause
%     end
    
    storePos = [];
    for s = 1:length(secs)
        tempID = str2num(sprintf('%02.0f%04.0f',w,secs(s)));
        tempPos = find(secID == tempID);
        if ~isempty(tempPos)
            storePos = [storePos tempPos];
        end
    end
    
    if length(storePos)>1
        sortPos = sort(storePos,'ascend');
        secOrder(sortPos) = secID(storePos);
    end
    
end
offOne = find(secOrder(2:end)-secOrder(1:end-1)~= 1);
showSwitches = [secOrder(offOne)' secOrder(offOne+1)'];

%%

for i = 1:length(secOrder)
    i
    targ = find(secID == secOrder(i));
    oldName = tifs(targ).name;
    newName = sprintf('z%05.f_%s',i,oldName);
   I(i,:) = imread([SPN 
    % copyfile([SPN oldName],[TPN newName]);
    
end
    
    


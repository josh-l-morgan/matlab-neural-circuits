
SPN = 'F:\S1-1Kx1Kcore_one\'
TPN = 'F:\S1-1Kx1Kcore_one_reorder\'
if ~exist(TPN,'dir'), mkdir(TPN),end

%nameGuide = {}

%% Parse nameGuide list of appropriate names
for i = 1:length(nameGuide);
    nam = nameGuide{i};
    zSec(i,:) = sscanf(nam,'z%d_W%d_s%d', [1, 3]);
    if ~isempty(regexp(nam,'dummy'))
        isDummy(i) = 1;
    else
        isdummy(i) = 0;
    end
end

%% Parse file names

oldNames = dir([SPN '*.tif']);
for i = 1:length(oldNames)
   nam = oldNames(i).name; 
   oldSec(i,:) = sscanf(nam,'W%d_s%d.tif',[1,2]);
end

%%

oldID = oldSec(:,1) * 1000 + oldSec(:,2);
newID = zSec(:,2) * 1000 + zSec(:,3);

isUsed = zeros(length(oldID),1);
for i = 1:length(newID)
   targ = find(oldID == newID(i));
   if isempty(targ)
       getImage(i) = 0;
   elseif length(targ)>1
       'too many targ found'
       break
   else
       getImage(i) = targ;
       isUsed(targ) = 1;
   end      
end
notUsed = sum(isUsed==0)


for i = 1:length(oldID)
   targ = find(newID == oldID(i));
   if isempty(targ)
       'missing targ'
       break
   elseif length(targ)>1
       'too many targ found'
       break
   else
       newOrder(i) = targ;
   end      
end



%%  write files
clear oldFile newFile
for i = 1:length(newOrder)
    oldFile{i,1} = oldNames(i).name;
    newFile{i,1} = sprintf('z%05.0f_W%02.0f_s%03.0f.tif',newOrder(i),oldSec(i,1),oldSec(i,2)); 
end

for i = 1:length(oldFile)
    i
    if ~exist([TPN newFile{i}],'file')
    copyfile([SPN oldFile{i}],[TPN newFile{i}])
    end 
end


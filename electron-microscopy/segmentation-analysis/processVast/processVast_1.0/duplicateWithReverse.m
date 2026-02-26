pause
SPN = GetMyDir;
dSPN = dir([SPN '*.tif'])
TPN =[SPN(1:end-1) 'withReverse\'];
mkdir(TPN)

%%
L = length(dSPN);
basename = dSPN(1).name(1:end-8);

for i = 1:L
    nam = dSPN(i).name;
    tifInd = str2num(nam(end-7:end-4));
    newname = sprintf('%s%04.0f.tif',basename,tifInd);
    copyfile([SPN nam],[TPN newname]);
    if (i <L) && ( i>1)
    newInd = L+ L-i;
    newername = sprintf('%s%04.0f.tif',basename,newInd);
    copyfile([SPN nam],[TPN newername]);

    end
end
    

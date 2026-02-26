

SPN = 'F:\samp10-24\W24\'
TPN = 'F:\samp10-24_one\'
mkdir(TPN)


dSPN = dir(SPN);
wafs = {};
for i = 1:length(dSPN)
    waf = dSPN(i).name;
    %if lower(waf(1)) == 'w'
        wafs{length(wafs)+1} = waf;
    %end
end

iNum = 0;
clear oldName newName
for w = 1:length(wafs)
    tWaf = dir([SPN wafs{w} '\*.tif']);
    for t = 1:length(tWaf)
        iNum = iNum+1;
        oldName{iNum} = sprintf('%s%s\\%s',SPN,wafs{w},tWaf(t).name);
        newName{iNum} = sprintf('%s%s_%s',TPN,wafs{w},tWaf(t).name);
    end
end


for i = 1:length(oldName)
    if ~exist(newName{i},'file')
        for r = 1:3
            pass = 1;
            try copyfile(oldName{i},newName{i})
            catch err
                pass = 0
            end
            if pass
                break
            end
        end
        
    end
    if ~mod(i,10)
        fileNum = length(dir([TPN '*.tif']));
        disp(sprintf('Copied %d of %d',fileNum,length(oldName)))
        pause(.01)
    end
end
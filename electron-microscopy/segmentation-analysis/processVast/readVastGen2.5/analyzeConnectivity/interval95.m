function[dataRange] = interval95(dat,intVal)

dat = dat(:);
if ~exist('intVal','var')
    intVal = 0.95;
end

step = (1-intVal)/2;

sortDat = sort(dat,'ascend');
datSize = length(dat);
lowRange = round(datSize*step)+1;
highRange = round(datSize - datSize*step)-1;

try
    dataRange = [sortDat(lowRange) sortDat(highRange)];
catch err
    err
    disp('bad interval')
end

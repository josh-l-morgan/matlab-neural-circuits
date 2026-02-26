
load('..\matFiles\us.mat')
%%
numSec = length(us.sec);
pathList = cell(numSec,16);
for s = 1:length(us.sec)
    for r = 1:4
        for c = 1:4
            pathList{s,(r-1)*4+c} = us.sec(s).paths{r,c};
        end
    end
end

xlswrite('..\xlsDir\pathList.xls','pathList')
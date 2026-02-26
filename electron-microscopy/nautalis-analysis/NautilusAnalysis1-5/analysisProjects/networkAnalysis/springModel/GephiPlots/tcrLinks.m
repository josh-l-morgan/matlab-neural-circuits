%function[csvText] = use2gephi(useList)


%% make node list
% 
% lineNum = 1
% colHead = 



mat2gephDir = 'D:\LGNs1\Analysis\connectivitySpreadsheets\'
nodeFileName = [mat2gephDir 'linkedTCRs.csv']
edgesFileName = [mat2gephDir 'links.csv']

useList = obI2cellList_seedInput(obI,[108 201])
conPref = seedPreferences([108 201],useList);



%%

lineNum = 1;
clear N
N{1,1} = 'ID';
N{1,2} = 'Label';
N{1,3} = 'Color';

for i = 1:length(useList.postList)
    lineNum = lineNum+1;
    N{lineNum,1} = num2str(useList.postList(i));
    N{lineNum,2} = 'tcr';
    N{lineNum,3} = '[1 0 0]';
end



delete(nodeFileName)
[nrows ncols] = size(N);
fid = fopen(nodeFileName, 'w');
for row=1:nrows
    lineStr = char();
    for col = 1:ncols
        lineStr = [lineStr N{row,col} ','];
    end
    lineStr = lineStr(1:end-1);
    fprintf(fid, '%s\n',lineStr);
end
fclose(fid);


%%



clear E
E{1,1} = 'Source';
E{1,2} = 'Target';
E{1,3} = 'Weight';
E{1,4} = 'Type';

con = useList.con;
linkNum = 1;
for x = 1:length(useList.postList)-1;
    for x2 = x+1:length(useList.postList);
        linkStrength =  sum(sqrt(con(:,x).* con(:,x2)));
        if linkStrength>0
        linkNum = linkNum+1;
        E{linkNum,1} = num2str(useList.postList(x));
        E{linkNum,2} = num2str(useList.postList(x2));
        E{linkNum,3} = num2str(linkStrength);
        E{linkNum,4} = 'Undirected';
        end
    end
end

delete(edgesFileName)
[nrows ncols] = size(E);
fid = fopen(edgesFileName, 'w');
for row=1:nrows
    lineStr = char();
    for col = 1:ncols
        lineStr = [lineStr E{row,col} ','];
    end
    lineStr = lineStr(1:end-1);
    fprintf(fid, '%s\n',lineStr);
end
fclose(fid);






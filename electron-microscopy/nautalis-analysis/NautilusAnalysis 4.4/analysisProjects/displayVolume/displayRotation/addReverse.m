function[] = addReverse(rotDir)





iDir = dir([rotDir '*.png']);
iNam = cat(1,{iDir.name});

newDir = [rotDir(1:end-1) 'reverse\'];
if~exist(newDir,'dir'), mkdir(newDir),end


for i  =1:length(iNam)
    nam = iNam{i};
    unds = find(nam == '_');
    iNum(i) = str2num(nam(unds(1)+1:unds(2)-1));
end
maxNum = max(iNum);

for i = 1:length(iNam)
    nam = iNam{i};
    unds = find(nam == '_');
    oldNum = str2num(nam(unds(1)+1:unds(2)-1));
    
    newNum = maxNum + maxNum - oldNum ;
    
    name1 = sprintf('%srot_%05.0f_.png',newDir,oldNum);
    name2 = sprintf('%srot_%05.0f_.png',newDir,newNum);
    oldName = [rotDir nam];
    copyfile(oldName,name1);
    copyfile(oldName,name2);
    
end



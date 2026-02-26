
FPN = GetMyDir;
ilist = dir([FPN '*.tif'])

qual = zeros(length(ilist),1);
parfor i = 1:length(ilist)
    disp(sprintf('checking %d of %d',i,length(ilist)))
    [tile maxId] = checkFileQualFull([FPN ilist(i).name]);
    qual(i)  = tile.quality;
    allQual(i).tile = tile;
    
end

% save('D:\LGNs1\Analysis\emVolume\a1kstackAllQual.mat','allQual')
%%
plot(qual)

lookRange = 8;
clear fixQual
for i = 1:length(allQual)
    fixQual(i) = allQual(i).tile.vertQual;
end
fixQual(isnan(fixQual)) = 0;



subplot(2,1,1)
plot(fixQual)
subplot(2,1,2)


for i = 1:length(qual)
    
    
    getWin= [i-lookRange:i-1 i+1:i+lookRange];
    getWin = getWin(getWin>=1);
    getWin =getWin(getWin<=length(qual));
    
    relqual(i) = fixQual(i)-median(fixQual(getWin));

end

plot(relqual)

stdRel = std(relqual);
upSTD = std(relqual(relqual>0))*3;
downSTD = std(relqual(relqual<0))*3;
higher = find(relqual>upSTD);
lower = find(relqual< - downSTD);

[sortRelQual sortIDX] = sort(relqual,'descend');

%plot(sortIDX)

sortIDX(1:10)
sortIDX(end-9:end)











%%


tot = 1701;

checked = 241;
checked = 700;

syn = 70;
syn = 300;

cells = 58;
cells = 200;


totSyn = round(tot * syn / checked);


checkNum = 1:totSyn;

reps = 1000;
clear medCount lowCount highCount
allCount = zeros(length(checkNum),reps);
for i = 1:length(checkNum)
   uCount = zeros(1,reps);
    for r = 1:reps
    mont = ceil(rand(totSyn,1) * checkNum(i));
    pick = randperm(totSyn,syn);
    allCell = mont(pick);
    uCell = unique(allCell);
    uCount(r) = length(uCell);
    end
    
    allCount(i,:) = uCount;
    
    sortCount = sort(uCount);
    medCount(i) = sortCount(round(reps/2));
    lowCount(i) = sortCount(round(reps*.05));
    highCount(i) = sortCount(round(reps*.95));
    
end


%%

plot(checkNum,medCount,'b')
hold on
plot(checkNum,lowCount,'r')
plot(checkNum,highCount,'r')
plot(checkNum,repmat(cells,[length(checkNum) 1]),'g')
hold off

%%

[y x ] = find(allCount == cells);
foundNum = sort(checkNum(y));
histFound = hist(foundNum,checkNum);
bar(checkNum,histFound)

lowFound = foundNum(round(length(foundNum)*.05))
midFound = foundNum(round(length(foundNum)*.5))
highFound = foundNum(round(length(foundNum)*.95))






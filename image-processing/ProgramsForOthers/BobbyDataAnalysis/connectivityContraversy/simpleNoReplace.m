opNum = 100;
axNum = 600;
synNum = 600;

reps = 1000;
histX = 0:1:2;

preIn = [];
for i = 1:axNum
    preIn = cat(1,preIn,ones(opNum,1)*i);
end

postIn = preIn * 0;
postIn(1:synNum) = 1;

allHist = zeros(reps,length(histX));
for r = 1:1000
  
    pos = rand(length(preIn),1);
    [a newPos] = sort(pos);
    newDend = postIn(newPos);
    
    con = preIn(newDend>0);
    
  sumCon = hist(con,1:1:axNum);
  histCon = hist(sumCon,histX);
  allHist(r,:) = histCon;     
    
end

meanHist = mean(allHist,1);
subplot(2,1,1)
bar(histX,meanHist)
subplot(2,1,2)
bar(meanHist(1:end-1)./meanHist(2:end))

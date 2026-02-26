opNum = 40;
axNum = 600;
synapseNumber = 180;

probSyn = synapseNumber/(opNum*axNum);

for r = 1:1000
  con = ( rand(axNum,opNum)) < probSyn  ;
  sumCon = sum(con,2);
  histCon = hist(sumCon,0:1:3);
  allHist(r,:) = histCon;  
    
    
end

meanHist = mean(allHist,1);
bar([0:1:3],meanHist)
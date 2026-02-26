


%%

%%Cell B
tc = 82
avSyn = 428
inNum = 17

%%CellA
tc = 68
avSyn = 212
inNum = 5


%%Cell E
tc = 19
avSyn = 24
inNum = 5

%%Type Mean
meanGiant = mean([62 62 26])
avSyn = mean([212 428  meanGiant ]); %seed cell average
avSyn = mean([212 428  62 62 26 ]); %seed cell average

avSyn = 50
%%

axBut = 125;

(tc * avSyn/inNum )/inNum + avSyn/inNum

divergence = axBut / (avSyn/inNum)
synSplit = avSyn/inNum
diverge  = axBut/synSplit
conVerge = avSyn/synSplit
synSplit = axBut/diverge
conVerge = avSyn/(axBut/diverge)
conVerge = avSyn/(axBut/(axBut/synSplit))

conVerge = avSyn/(axBut/(axBut/(avSyn/inNum)))



%%Bouton per axon estimate baseed on particular in num
testNum = [1:30];
butPerAx = tc*avSyn/inNum/inNum + avSyn/inNum
testButPerAx = tc*avSyn./testNum./testNum + avSyn./testNum
bestNum = testNum(max(find(testButPerAx>axBut)))


%{
butPerAx/aveSyn = tc * 1/inNum * 1/inNum + 1/inNum;
butPerAx * inNum * 1/inNum = tc * inNum + 1;
butPerAx = tc*inNum+1;
%}

inNum = (axBut - 1)/tc 


%%Convergence based on reported bouton number
converge = tc*avSyn/axBut

converge = sqrt(tc*avSyn/axBut)

%%


avPreSyn = 125;
avPostSyn = 150;
avSubnetSize = 52;







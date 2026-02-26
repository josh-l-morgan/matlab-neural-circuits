function[result] = scatterMonte(resReal, resMonte,resName)


monteMin = min(resMonte);
monteMax = max(resMonte);
monteRange = range(resMonte);

datMin = min(min(resReal),monteMin);
datMax = max(max(resReal),monteMax);
datRange = datMax-datMin;

[histMonte xout] = hist(resMonte,[monteMin:monteRange/10:monteMax]);
maxHist = max(histMonte);


bar(xout, histMonte,'barWidth',1)
hold on
scatter(resReal,1,'r')
xlim([datMin-datRange/10 datMax+datRange/10 ])
hold off
intVal = interval95(resMonte);
Pgreater = sum((resMonte>=resReal)/length(resMonte));

disText = sprintf('%s -> real = [%.2f] \n monte 95range = [%.2f - %.2f]\n P = %.7f',...
    resName,resReal,intVal(1),intVal(2),Pgreater);


if resReal > mean(resMonte)
    text(datMax-datRange/3,maxHist-maxHist/5, disText)
else
    text(datMin,maxHist-maxHist/5, disText)
end


disp(sprintf('%s -> real = [%.2f], monte 95range = [%.2f - %.2f], P = %.7f',...
    resName,resReal,intVal(1),intVal(2),Pgreater))

sortMonte = sort(resMonte,'ascend');

result.name = resName;
result.real = resReal;
result.monte95range = intVal;
result.monteReps = length(resMonte);
result.Pgreater = Pgreater;


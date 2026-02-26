function[h] = barScat2(varargin)

Yin = varargin;
binWidth = 3;

Xax = 1: length(Yin);
theMax = 0;
theMin = 0;
for i = 1:length(Yin);
    Y = Yin{i};
    meanY(i) = median(Y);
    Ns(i) = length(Y);
    E(i) = std(Y)/sqrt(Ns(i));
    theMax = max(theMax,max(Y));
    theMin = min(theMin,min(Y));
end

bar(meanY,'w')
hold on 

h = errorbar(Xax,meanY,E,'k','LineWidth',2,'MarkerSize',10);
errorbar_tick(h,3 * max(Xax));
for i = 1:length(Yin)
    scatY = Yin{i};
    scatX = ones(length(Yin{i}),1)*i;
    
    %% bump
    xBin = min(scatY):(std(scatY)/binWidth):(max(scatY)+std(scatY)/binWidth);
    xHist = histc(scatY,xBin);
    maxH = max(xHist);
    for h = 1 : length(xBin)-1
       tobump = find((scatY >= xBin(h)) & (scatY < xBin(h+1)));
       numbump = length(tobump);
       bumper = (1:numbump)* .4/maxH;
       bumper = bumper - mean(bumper);
       scatX(tobump) = scatX(tobump) + bumper';
    end

    
    scatter(scatX,scatY,'k','LineWidth',1);
    text(i,theMax+theMax/5,num2str(Ns(i)))
end

ylim([theMin - theMin/4 theMax + theMax/4])

hold off
% 
% ys=ylim;
% ylim([min(0,ys(1)) ys(2)]);
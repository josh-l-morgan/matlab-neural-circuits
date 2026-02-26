function[h] = barScat5(varargin)

Yin = varargin;
binWidth = 1;

Xax = 1: length(Yin);
theMax = 0;
theMin = 0;
for i = 1:length(Yin);
    Y = Yin{i};
    meanY(i) = mean(Y);
    Ns(i) = length(Y);
    E(i) = std(Y)/sqrt(Ns(i));
    theMax = max(theMax,max(Y));
    theMin = min(theMin,min(Y));
end

bar(meanY,'FaceColor',[.8 .8 .8],'BarWidth',.9)
hold on 
% 
% h = errorbar(Xax,meanY,E,'k','LineWidth',2,'MarkerSize',10);
% errorbar_tick(h,3 * max(Xax));
for i = 1:length(Yin)
    scatY = Yin{i};
    scatX = ones(length(Yin{i}),1)*i;
    
    %% bump
    binStep = .001;%8;%std(scatY)/binWidth;
    binStep = (max(scatY) - min(scatY))/100;
    if ~binStep, binStep = 0.001; end
    xBin = min(scatY): binStep :max(scatY) + binStep;
    
    xHist = histc(scatY,xBin);
    maxH = max(xHist);
    maxH = max(4,maxH);
    for h = 1 : length(xBin)-1
       tobump = find((scatY >= xBin(h)) & (scatY < xBin(h+1)));
       numbump = length(tobump);
       bumper = (1:numbump)* .7/maxH;
       bumper = bumper - mean(bumper);
       scatX(tobump) = scatX(tobump) + bumper';
    end

    
    scatter(scatX,scatY,'k','LineWidth',1);
    xlim([.45 2.55])
%     text(i,theMax+theMax/5,num2str(Ns(i)))
end

ylim([theMin - theMin/4 theMax + theMax/4])
P = ranksum(Yin{1},Yin{2});
% if P > 0.001
% xlabel(sprintf('P = %.3f',P));
% else
%     xlabel('P < 0.001')
% end

hold off
% 
% ys=ylim;
% ylim([min(0,ys(1)) ys(2)]);
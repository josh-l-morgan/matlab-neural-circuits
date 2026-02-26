function[h] = barScat(varargin)

Xax = 1: length(varargin);
theMax = 0;
for i = 1:length(varargin);
    Y = varargin{i};
    meanY(i) = mean(Y);
    Ns(i) = length(Y);
    E(i) = std(Y)/sqrt(Ns(i));
    theMax = max(theMax,max(Y));
end

bar(meanY,'w')
hold on 

h = errorbar(Xax,meanY,E,'k','LineWidth',2,'MarkerSize',10);
errorbar_tick(h,3 * max(Xax));
for i = 1:length(varargin)
    scatY = varargin{i};
    scatX = ones(length(varargin{i}),1)*i;
    %hScat = hist(scatY,0:max(scatY));
    
    
    
    scatter(scatX,scatY,'k','LineWidth',2);
    text(i,theMax+theMax/5,num2str(Ns(i)))
end

ylim([0 theMax + theMax/4])

hold off
% 
% ys=ylim;
% ylim([min(0,ys(1)) ys(2)]);
function[] = plotBin(X,Y,B,C)

if isempty(X) | isempty(Y)
    'no data to plot'
    plot(0,0,'w')
else

    if nargin ==2
        C = 'b';
        B = 1;
    elseif nargin == 3
        C = 'b';
    end


    X = X(:); Y = Y(:);
    step = (fix(max(X))- fix(min(X)))/30;
    Xax = fix(min(X)) : B/5 : fix(max(X))+1 + B/10;
    Xax = B/2:B:1-B/2
    Yold = [];P = [];
    for i = 1:length(Xax);
        
        use = (X >= Xax(i)-B/2) & (X < Xax(i)+B/2);
        meanY(i,1) = mean(Y(use));
        E(i) = std(Y(use))/sqrt(sum(use));
        if ~isempty(Yold)
            P(i) = ranksum(Yold,Y(use));
        end
        Yold = Y(use);
    end

    %errorbar(Xax,meanY,E,C,'LineWidth',1,'MarkerSize',440)
    bar(Xax,meanY,'BarWidth',1,'FaceColor',[.8 .8 .8])
    hold on
    scatter(X,Y,C,'LineWidth',1)
    legend(num2str(P))
    hold off
    %
    % ys=ylim;
    % ylim([min(0,ys(1)) ys(2)]);
end
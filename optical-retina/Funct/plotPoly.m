function[Slope Std] = plotPoly(X,Y,C)

if isempty(X) | isempty(Y)
    'no data to plot'
    plot(0,0,'w')
else

    if nargin ==2
        C = 'b';
    end


    X = X(:); Y = Y(:);
    Xn = (X-mean(X))/std(X); Yn = (Y-mean(Y))/std(Y);
    [p, S, mu] = polyfit(Xn,Yn,1);
    f = polyval(p,Xn) * std(Y) + mean(Y);
    p
    table = [X Y f Y-f];
    SE = std(Y-f)/sqrt(length(Y));
    Std = std(Y-f);

    
    idXmin = find(X == min(X),1);
    idXmax = find(X == max(X),1);
    
    Slope(2) = (f(idXmax)-f(idXmin))/(X(idXmax)-X(idXmin));
    Slope(1) = polyval(p,0) * std(Y) + mean(Y);
    
    plot(X,f,C,'LineWidth',1);
    hold on
    scatter(X,Y,C,'LineWidth',2)
    hold off
    %
    % ys=ylim;
    % ylim([min(0,ys(1)) ys(2)]);
%     [R P] = corrcoef(Xn, Yn);
%     R = R(1,2); P = P(1,2);
    
end
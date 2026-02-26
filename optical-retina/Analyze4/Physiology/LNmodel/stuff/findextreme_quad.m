function ind_quad = findextreme_quad(input)
%this function returns the linear index of the extremum of a vector or
%matrix

if max(input) > abs(min(input))
    ind = find(input == max(input));
    y = input((ind-1):(ind+1));
    x = -1:1;
    p = polyfit(x,y,2);
    xfit = linspace(-1,1,20);
    yfit = polyval(p, xfit);
    ind_quad = ind + xfit(yfit == max(yfit));

elseif max(input) < abs(min(input))
    ind = find(input == min(input));
    y = input((ind-1):(ind+1));
    x = -1:1;
    p = polyfit(x,y,2);
    xfit = linspace(-1,1,20);
    yfit = polyval(p, xfit);
    ind_quad = ind + xfit(yfit == min(yfit));
else
    disp ('NO UNIQUE EXRTEMUM')
end

if length(ind) > 1
    disp ('NO UNIQUE EXTREMUM')
else
end
function plotxyerrorD(x,y,xe,ye,marker,line)

%         plotxyerrorD(x,y,xe,ye,marker,line)
%         x,y =     [n x 1] vectors of data values.
%         xe =      [n x 1] corresponding vector of horizontal deviations
%         ye =      [n x 1] corresponding vector of vertical deviations
%         marker = string defining markersymbol and color (e.g. 'rs')
%         line = string defining line color (e.g. 'r')

n = length(x);
plot(x,y,marker)
hold on
for i=1:n
       plot([x(i)-xe(i) x(i)+xe(i)],[y(i) y(i)],line)
       plot([x(i) x(i)],[y(i)-ye(i) y(i)+ye(i)],line)
end



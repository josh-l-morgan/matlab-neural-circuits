function drawpoly(x,y,color)
%draws a polygon using the (x,y) coordinate vectors
%by Daniel Berger for MIT/BCS Seung, April 21 2009

xs=size(x);
ys=size(y);
if (size(xs,2)==2) && (xs(1)==ys(1)) && (xs(2)==ys(2))
  
  if xs(1)<xs(2) %make column vectors
    x=x'; y=y';
  end;
  
  if size(x,2)==1 %only if column vector
    plot(x,y,color,'linewidth',2);
    if (ishold==0)
      hold on;
      plot([x(size(x,1)) x(1)],[y(size(y,1)) y(1)],color,'linewidth',2);
      hold off;
    else
      plot([x(size(x,1)) x(1)],[y(size(y,1)) y(1)],color,'linewidth',2);  
    end;
  end;
end;
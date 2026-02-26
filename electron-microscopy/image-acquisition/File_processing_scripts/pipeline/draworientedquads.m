function draworientedquads(sgridx,sgridy,radang,xedgelength,yedgelength,color)

cornerx=zeros(5,1);
cornery=zeros(5,1);
% [cornerx(1),cornery(1)]=rotatevec(-1,-1,radang);
% [cornerx(2),cornery(2)]=rotatevec(-1,1,radang);
% [cornerx(3),cornery(3)]=rotatevec(1,1,radang);
% [cornerx(4),cornery(4)]=rotatevec(1,-1,radang);
% [cornerx(5),cornery(5)]=rotatevec(-1,-1,radang);
% cornerx=cornerx*xedgelength/2;
% cornery=cornery*yedgelength/2;

xh=xedgelength/2; yh=yedgelength/2;
[cornerx(1),cornery(1)]=rotatevec(-xh,-yh,radang);
[cornerx(2),cornery(2)]=rotatevec(-xh,yh,radang);
[cornerx(3),cornery(3)]=rotatevec(xh,yh,radang);
[cornerx(4),cornery(4)]=rotatevec(xh,-yh,radang);
[cornerx(5),cornery(5)]=rotatevec(-xh,-yh,radang);

for y=1:1:size(sgridx,1)
  for x=1:1:size(sgridx,2)
    xp=sgridx(y,x);
    yp=sgridy(y,x);
    plot(cornerx(1)+xp,cornery(1)+yp,'k.');
    plot(cornerx+xp,cornery+yp,color);
  end;
end;
function draworientedsquares(sgridx,sgridy,radang,edgelength,color)

cornerx=zeros(5,1);
cornery=zeros(5,1);
[cornerx(1),cornery(1)]=rotatevec(-1,-1,radang);
[cornerx(2),cornery(2)]=rotatevec(-1,1,radang);
[cornerx(3),cornery(3)]=rotatevec(1,1,radang);
[cornerx(4),cornery(4)]=rotatevec(1,-1,radang);
[cornerx(5),cornery(5)]=rotatevec(-1,-1,radang);
cornerx=cornerx*edgelength/2;
cornery=cornery*edgelength/2;

for y=1:1:size(sgridx,1)
  for x=1:1:size(sgridx,2)
    xp=sgridx(y,x);
    yp=sgridy(y,x);
    plot(cornerx(1)+xp,cornery(1)+yp,'k.');
    plot(cornerx+xp,cornery+yp,color);
  end;
end;
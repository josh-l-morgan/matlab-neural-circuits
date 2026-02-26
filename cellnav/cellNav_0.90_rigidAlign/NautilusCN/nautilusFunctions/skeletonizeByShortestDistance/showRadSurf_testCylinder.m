




pos = [ 0 0 0; 1 -1 -1; .5 .5 .5 ; .5 2 1];
edge = [1 2; 2 3; 3 4]
rad = [.1 .3 .2 .1];


div = pi/10;
d = div:div:(2 * pi);
x = cos(d)';
y = sin(d)';
z = x * 0';

cSub = [x y z x*0 ]

%%
clf
hold on

daspect([1,1,1])
view(3); axis tight
set(gca,'color',[0 0 0])
set(gcf,'color',[0 0 0])
%set(gca,'off')
l = lightangle(145,45) ;


for i = 1:1;%size(edge,1)
    
    i = 1
    p1 =  pos(edge(i,1),:);
    p2 = pos(edge(i,2),:);
    dif =   p2 - p1;
    ar1= atan2(dif(1),dif(2))-pi/2;
    ax = [sin(ar1) cos(ar1) 0];
    dist = sqrt(sum(dif.^2));
    ar2 = asin(dif(3)/dist)+pi/2;
%     axDist = sqrt(dif(2)^2+dif(1)^2);
%     ar2 = atan2(axDist,dif(3));
%     
%     r1 = atan2(dif(1),dif(2));
%     r2 = atan2(dif(3),dif(1));
%     
%     mz = makehgtform('zrotate',zr(a));
%     mx = makehgtform('xrotate',xr(a));
    mt = makehgtform('axisrotate',ax,ar2);
    
    %mt = mz * mx;
    
    
%     E = difs(a,:);
%     %Direction Cosines (rotation matrix) construction
%     Ry=[1        0        0;...
%         0        cos(E(1))  -sin(E(1));...
%         0        sin(E(1))  cos(E(1))]; %X-Axis rotation
%     Rx=[cos(E(2))  0        sin(E(2));...
%         0        1        0;...
%         -sin(E(2)) 0        cos(E(2))]; %Y-axis rotation
%     Rz=[cos(E(3))  -sin(E(3)) 0;...
%         sin(E(3))  cos(E(3))  0;...
%         0        0        1]; %Z-axis rotation
%     R=Rx*Ry*Rz; %Rotation matrix
%     
%     % mt = R;
    
    rSub1 = cSub * mt;
    rSub1 = rSub1 * rad(edge(i,1));
    rSub1 = rSub1(:,1:3) + repmat(p1,[size(rSub1,1) 1]);
    rSub2 = cSub * mt;
    rSub2 = rSub2(:,1:3)  * rad(edge(i,2));
    rSub2 = rSub2 + repmat(p2,[size(rSub1,1) 1]);
    
    plot3(rSub1(:,1),rSub1(:,2),rSub1(:,3),'r')
    plot3(rSub2(:,1),rSub2(:,2),rSub2(:,3),'g')
    plot3([0 ax(1)],[0 ax(2)],[0 ax(3)],'w')
    plot3([0 dif(1)],[0 dif(2)],[0 dif(3)],'y')

    pause(1)
    
end

%end

%%
























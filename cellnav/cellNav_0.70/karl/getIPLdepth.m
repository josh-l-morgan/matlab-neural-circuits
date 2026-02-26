%get the depth of a point in the IPL from the fitted planes
function [zGCL,zINL,IPLdepth] = getIPLdepth(z,x,y,GCLplane,INLplane)
zGCL=(-double(x)*GCLplane.Parameters(2)-double(y)*GCLplane.Parameters(3)-GCLplane.Parameters(4))/GCLplane.Parameters(1);
zINL=(-double(x)*INLplane.Parameters(2)-double(y)*INLplane.Parameters(3)-INLplane.Parameters(4))/INLplane.Parameters(1);
absBool=0;
if absBool==1
    IPLdepth=abs(double(z)-zINL)./abs(zINL-zGCL);
else
    IPLdepth=(zINL-double(z))./(zINL-zGCL);
end
if 1 %eyewireFit==1
    m=0.9075;
    b=0.0052;
    IPLdepth=IPLdepth.*m+b;
end
end
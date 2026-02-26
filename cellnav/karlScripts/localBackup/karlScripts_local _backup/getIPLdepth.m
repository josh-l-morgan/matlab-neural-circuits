function [zGCL,zINL,IPLdepth] = getIPLdepth(z,x,y,GCLplane,INLplane)
%if using locations from tis.mat, use location([3 1 2])
%if using FV locations, use location([3 2 1]);
if isempty(GCLplane)
    GCLplane=struct();
    GCLplane.Parameters=[-1.0000   -0.0525    0.0262   10.3135];
end
if isempty(INLplane)
    INLplane=struct();
    INLplane.Parameters=[-1.0000   -0.0768    0.0327   43.0766];
end
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
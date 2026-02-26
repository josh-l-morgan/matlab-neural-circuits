function outputZs=getZdistFromUm(inputCoords)

function [zGCL,zINL,IPLdepth] = getIPLdepth(z,x,y,GCLplane,INLplane)
zGCL=(-x*GCLplane.Parameters(2)-y*GCLplane.Parameters(3)-GCLplane.Parameters(4))/GCLplane.Parameters(1);
zINL=(-x*INLplane.Parameters(2)-y*INLplane.Parameters(3)-INLplane.Parameters(4))/INLplane.Parameters(1);

IPLdepth=abs(z-zINL)/abs(zINL-zGCL);
if 1 %eyewireFit==1
    m=0.9075;
    b=0.0052;
    IPLdepth=IPLdepth*m+b;
end
end



end
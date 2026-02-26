function[newI] = buffI(I,buf)
%%adds or subtracts buf number of voxesl to the front and back of each
%%dimension of a three dimentional stack

[ys xs zs] = size(I);
if buf >0
    newI = zeros(ys+2*buf, xs+2*buf,zs+2*buf);
    newI(buf+1:buf+ys,buf+1:buf+xs,buf+1:buf+zs) = I;
else
    newI = I(1-buf:ys+buf,1-buf:xs+buf,1-buf:zs+buf);
end

function[I2] = pad(I,buf)

I2 = zeros(size(I,1)+buf * 2, size(I,2) + buf * 2,class(I));
I2(buf+1:end-buf,buf+1:end-buf) = I;
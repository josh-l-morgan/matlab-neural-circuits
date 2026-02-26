function[I2] = unpad(I,buf)

I2 = I(buf + 1:end-buf,buf +1:end-buf);
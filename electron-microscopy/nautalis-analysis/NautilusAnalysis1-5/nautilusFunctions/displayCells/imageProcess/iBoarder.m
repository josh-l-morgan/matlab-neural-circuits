function[I] = iBoarder(I,buf,bufVal)



I(1:buf,:,:) = bufVal;
I(end-buf:end,:,:) = bufVal;
I(:,1:buf,:) = bufVal;
I(:,end-buf:end,:,:) = bufVal;





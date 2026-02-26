 
 MIJ.start('C:\Program Files\Fiji.app')
 
 fsize = 200;
 gKern = gaus3d([10 10 1],3);
 I1 = rand(fsize);
 I2 = circshift(I1,10);
 I1 = fastCon(I1,gKern);
 I2 = fastCon(I2,gKern);
 
[X,Y,R] = DetermineAlignmentUsingSIFT(uint8(I1),uint8(I2))
%%
topSamp = top;
image(topSamp)
ISamp = I(1000:3000,:);
image(ISamp)
[X,Y,R] = DetermineAlignmentUsingSIFT(ISamp,topSamp)
[X,Y,R] = mijSIFT(ISamp,topSamp)


%%

image(top(1:fsize,1:fsize))



KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;

for k = 1:size(Kdir,1)
    TPN = [KPN Kdir(k).name ]
    I = imread(TPN);
    I(I>0)=255;    
   
    Name=[KPN Kdir(k).name];
    imwrite(I,Name,'Compression','none')

 
end
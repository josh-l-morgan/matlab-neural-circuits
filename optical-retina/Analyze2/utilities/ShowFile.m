clear all

KPN=GetMyDir
Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;

for k = 1:size(Kdir,1)
    TPN = [KPN '\' Kdir(k).name '\']; 
    if exist([TPN 'find\SG.mat'])
   
        load([TPN 'Use.mat'])
        Use
        
    end
end
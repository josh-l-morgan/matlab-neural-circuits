clear all

KPN=GetMyDir
Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;



load('cmap.mat')
set(gcf,'Colormap',cmap)

for k = 1:1%size(Kdir,1)
    clear Branch 
    TPN = [KPN '\' Kdir(k).name '\']
   if exist([TPN 'Branch.mat'])
       load([TPN 'Branch.mat'])
       
       
   end
end
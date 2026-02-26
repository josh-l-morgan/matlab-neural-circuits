%%Specifically finds overlap with listed ablated territory

clear all
KPN=GetMyDir
TPN=[KPN 'Ter\']
Idir=dir(TPN); Idir=Idir(3:size(Idir));
colormap gray(255)

Names={};
for k = 1 : size(Idir,1)
Names(k,1)={Idir(k).name};    
Ters(:,:,k)=imread([TPN Idir(k).name]);
%'Raw Image'
end
Ters=Ters>0;


slash=find(KPN=='\');
GPN=KPN(1:slash(size(slash,2)-1));
GPNf=[GPN 'PreUn\Ter'];
GPNd=dir(GPNf);
GPNd=GPNd(3:size(GPNd,1));
for i = 1: size(GPNd,1)
    UnablatedLocation=[GPNf '\' GPNd(i).name];
    UnabRegi=imread(UnablatedLocation);
    UnabRegi=UnabRegi>0;
    UnabRegi=max(UnabRegi,[],3);
    UnabReg(:,:,i)=UnabRegi;
end
image(sum(UnabReg,3)*50)


%%
SumTers=sum(Ters,3);
image(SumTers*50+sum(UnabReg,3)*30),pause(.01)

Dat(1,1)={'Area'};
Dat(1,2)={'Overlap'};
Dat(1,3)={'Private'};


for i = 1:size(Ters,3)
   Area(i)=sum(sum(Ters(:,:,i))) ;
   for t = 1 : size(Ters,3) %% run against neighbors
        Nieghbor(i,t)=(sum(sum(Ters(:,:,i)& UnabReg(:,:,t))))/Area(i);
   end    
   Overlap(i)=sum(sum(SumTers(Ters(:,:,i)>0)))-Area(i);
   Private(i)=Area(i)-Overlap(i); 
   
   Dat(i+1,1)={Area(i)};
   Dat(i+1,2)={Overlap(i)};
   Dat(i+1,3)={Private(i)};
   
end
Dat(2:size(Names,1)+1,4)=Names;
Dat
xlswrite([KPN 'FirstNeighborDat.xls'],Nieghbor)












     
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
GPNf=[GPN 'Combined\Ter'];
GPNd=dir(GPNf);
AblatedLocation=[GPNf '\' GPNd(3).name];
AbReg=imread(AblatedLocation);
AbReg=AbReg>0;
AbReg=max(AbReg,[],3);

%%
SumTers=sum(Ters,3);
image(SumTers*50+AbReg*100),pause(.01)

Dat(1,1)={'Area'};
Dat(1,2)={'Overlap'};
Dat(1,3)={'Private'};

for i = 1:size(Ters,3)
   Area(i)=sum(sum(Ters(:,:,i))) ;
    image(AbReg * 100 + Ters(:,:,i)*200)
    pause(.1)
   
   Overlap(i)=sum(sum(AbReg(Ters(:,:,i)>0)));
   Private(i)=Area(i)-Overlap(i); 
   
   Dat(i+1,1)={Area(i)};
   Dat(i+1,2)={Overlap(i)};
   Dat(i+1,3)={Private(i)};
   
end
Dat(2:size(Names,1)+1,4)=Names;
Dat
xlswrite([KPN 'AbDat.xls'],Dat)













     
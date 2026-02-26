
clear all

KPN=GetMyDir
Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;

for k = 1:1%size(Kdir,1)
    TPN = [KPN '\' Kdir(k).name '\']; 
   if exist([TPN 'mask'])
       
       
       %% Read Mask
       'Read Mask'
       IMdir=dir([TPN 'mask']); IMdir=IMdir(3:size(IMdir,1));
       IM=imread([TPN 'mask\' IMdir(1).name]);
       [ys xs zs] = size(IM);
       IM=zeros(ys,xs,size(IMdir,1),'uint8');
       for i = 1:size(IMdir,1)
          IM(:,:,i)=imread([TPN 'mask\' IMdir(i).name]); 
       end
       
       %% Find Volume
       'Find Volume'
       IM=IM==1;
       [y x z]=find3(IM);
       Msize=fix([ys xs]*.103)+1;
       %%Scale and correct Volume index
       y=round(y*.103);
       x=round(x*.103);
       y(y<1)=1;
       y(y>ys)=ys;
       x(x<1)=1;
       x(x>xs)=xs;
       
       imVol=zeros(Msize(1),Msize(2));
       for i = 1: size(y,1)
          imVol(y(i),x(i))=imVol(y(i),x(i))+yxum*yxum*zum;
       end
       image(imVol*300), pause(.3)
       clear y x z
       
       %% Find Surface
       'Find Surface'
       IM=bwperim(IM,18);
       [y x z]=find3(IM);
       Msize=fix([ys xs]*.103)+1;
       %%Scale and correct Volume index
       y=round(y*.103);
       x=round(x*.103);
       y(y<1)=1;
       y(y>ys)=ys;
       x(x<1)=1;
       x(x>xs)=xs; 
       
       imSurf=zeros(Msize(1),Msize(2));
       for i = 1: size(y,1)
          imSurf(y(i),x(i))=imSurf(y(i),x(i))+.03;
       end
       image(imSurf*300), pause(.3)
       clear y x z
       
    
    
   end
end
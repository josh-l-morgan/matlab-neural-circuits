
%clear all

KPN=GetMyDir
Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;
colormap(jet(256))
cmap=colormap;
cmap(1,:)=0;
colormap(cmap)



for k = 1:size(Kdir,1)
    clear CA4   
    TPN = [KPN '\' Kdir(k).name '\']; 
     clear CA BiDotDist DotShow DotFilt DendFilt DotMap DendMap DPos 
     clear Cell Red Green Blue TerFill TerFilt TerPerim TerRay Territory
    if exist([TPN 'CA.mat'])
      load([TPN 'CA.mat'])
      
      DotMap=CA.Arbor(1).DotMap;
      DendMap=CA.Arbor(1).DendMap;
      DotDist=CA.Arbor(1).DotDist;
      DendDist=CA.Arbor(1).DendDist;
      Territory =CA.Arbor(1).Territory;
      
      [ys xs]=size(DotMap);
      
      DotDend=DendMap;
      DotDend(DotDend>0)=245;
      DotDend(DotMap>0)=100;
      image(DotDend)
            
      CA4=zeros(ys,844);
      CA4(1:ys,423:423+xs-1)=DotDist*1500+Territory*2;
      CA4(1:ys,212:212+xs-1)=DendDist*300+Territory*2;
      CA4(1:ys,1:xs)=DotDend;
      DDi=DotDist./DendDist*300+Territory*2;
      DDi(DendDist==0)=Territory(DendDist==0)*2;
      CA4(1:ys,634:634+xs-1)=DDi;
      image(CA4),pause(.3)
      
      
      Perimeter=bwperim(Territory);
      
      
      if exist([TPN 'Cell.mat'])
            load([TPN 'Cell.mat'])
            cName=['P' Cell.Age '_' Cell.Type '_' Kdir(k).name '.tif']
            Name=[KPN 'Images\C4B\' cName];
            Name2=[KPN 'Images\Territories\' cName];
        else
            Name=[KPN 'Images\C4B\' Kdir(k).name '.tif'];
            Name2=[KPN 'Images\Territories\' Kdir(k).name '.tif'];
        end
        if isempty(find(Name=='?'));
            imwrite(CA4,cmap,Name,'Compression','none')
            imwrite(Perimeter,Name2,'Compression','none')
        end
      
        
        
        
      
      
          
    end
end



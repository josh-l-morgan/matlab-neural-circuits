
%clear all
colormap(jet(256))
cmap=colormap;
cmap(1,:)=0;
colormap(cmap)

load('UseCells.mat')
KPN = '\\128.208.64.36\wonglab\Josh\Analyzed\Output2\'

for k = 1:size(UseCells,1)
    clear CA4   
    TPN = char(UseCells(k)); 
     clear CA BiDotDist DotShow DotFilt DendFilt DotMap DendMap DPos 
     clear Cell Red Green Blue TerFill TerFilt TerPerim TerRay Territory
    if exist([TPN 'CA.mat'])
      load([TPN 'CASc.mat'])
      
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
            
      CA4=zeros(ys,961);
      CA4(1:ys,481:481+xs-1)=DotDist*1500+Territory*2;
      CA4(1:ys,241:241+xs-1)=DendDist*300+Territory*2;
      CA4(1:ys,1:xs)=DotDend;
      DDi=DotDist./DendDist*300+Territory*2;
      DDi(DendDist==0)=Territory(DendDist==0)*2;
      CA4(1:ys,721:721+xs-1)=DDi;
      image(CA4),pause(.3)
      
      
      Perimeter=bwperim(Territory);
      
      
      if exist([TPN 'Cell.mat'])
            load([TPN 'Cell.mat'])
            cName=['P' Cell.Age '_' Cell.Type '_' Cell.Name '.tif']
            Name=[KPN 'Images\C4B\' cName];
            Name2=[KPN 'Images\Territories\' cName];
        else
            Name=[KPN 'Images\C4B\' Kdir(k).name '.tif'];
            Name2=[KPN 'Images\Territories\' Cell.Name '.tif'];
        end
        if isempty(find(Name=='?'));
            imwrite(CA4,cmap,Name,'Compression','none')
            imwrite(Perimeter,Name2,'Compression','none')
        end
      
        
        
        
      
      
          
    end
end



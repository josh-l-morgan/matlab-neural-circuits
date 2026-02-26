function CAsampleCollectN2007(TPN)

clf
colormap(jet(256))
cmap=colormap;
cmap(1,:)=0;
colormap(cmap)

% load('UseCells.mat')
% KPN = '\\128.208.64.36\wonglab\Josh\Analyzed\Output2\'

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

    CA4=zeros(ys,961);
    CA4(1:ys,481:481+xs-1)=DotDist*1500+Territory*2;
    CA4(1:ys,241:241+xs-1)=DendDist*300+Territory*2;
    CA4(1:ys,1:xs)=DotDend;
    DDi=DotDist./DendDist*300+Territory*2;
    DDi(DendDist==0)=Territory(DendDist==0)*2;
    CA4(1:ys,721:721+xs-1)=DDi;
    image(CA4),pause(.3)


    Perimeter=bwperim(Territory);


    Name=[TPN 'Images\C4B\Bild'  '.tif'];
    Name2=[TPN 'Images\Territories\Bild' '.tif'];
    
    if isempty(find(Name=='?'));
        imwrite(CA4,cmap,Name,'Compression','none')
        imwrite(Perimeter,Name2,'Compression','none')
    end







end




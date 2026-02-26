%%Get Numbers for Inner Half, Outer Half, Inner Third, Outer Third, 
%%inner half - CB, Outer Half / Inner Half - CB

clear all


load('UseCells.mat')

II=.1; %inner boarder

for k = 1:1%size(UseCells,1)
    k
    TPN = char(UseCells(k));
    clear Grad
    
    %% Run as vectors
    load([TPN 'Use.mat'])
    load([TPN 'data\AllSegCut.mat'])
    Mids=mean(AllSegCut,3);
    Length=sqrt((AllSegCut(:,1,1)-AllSegCut(:,1,2)).^2 ...
           + (AllSegCut(:,2,1)-AllSegCut(:,2,2)).^2 ...
           + (AllSegCut(:,3,1)-AllSegCut(:,3,2)).^2);
    
    %Find eccentricities
    EccDPos=sqrt((Use.DPos(:,1)-Use.Cent(1)).^2 + (Use.DPos(:,2)-Use.Cent(2)).^2);
    EccMids=sqrt((Mids(:,1)-Use.Cent(1)).^2 + (Mids(:,2)-Use.Cent(2)).^2);
    
    %%find edges
    SortMids=sort(EccMids);
    Outer=SortMids(fix(size(SortMids,1)*.98));
    Grad.Outer=Outer;
    
    InnerLength=sum(Length((EccMids>10) & (EccMids<(Outer/2))));
    OuterLength=sum(Length((EccMids>(Outer/2))));
    InnerDots=sum((EccDPos>10) & (EccDPos<(Outer/2)));
    OuterDots=sum((EccDPos>(Outer/2)));  
    Grad.Vect=(InnerDots/InnerLength)/(OuterDots/OuterLength);
    GradV(k)=Grad.Vect;
    if GradV(k)<.1,pause,end
    clear AllSegCut Mids Length
    
    %% Run from areas
    load([TPN 'CA.mat'])
    if size(CA.Arbor,2)>2
        Area=CA.Arbor(2).Territory+CA.Arbor(3).Territory;
        DotsA=CA.Arbor(2).DotMap + CA.Arbor(3).DotMap;
        DendA=CA.Arbor(2).DendMap + CA.Arbor(3).DendMap;
        
    else
        Area=CA.Arbor(1).Territory;
        DotsA=CA.Arbor(1).DotMap;
        DendA=CA.Arbor(1).DendMap;
    end
    
    DistMap=zeros(size(Area));
    for y = 1:size(Area,1)
        for x = 1: size(Area,2)
            DistMap(y,x)=sqrt((y-Use.Cent(1))^2+(x-Use.Cent(2))^2);
        end
    end
    
    image(((DistMap>10) & (DistMap<(Outer)))*50+Area*50),pause(.01)

    MaskA=((DistMap>10) & (DistMap<Outer))
    Operim=bwperim(MaskA);
    image(Operim*1000)
    ShowOuter=uint8((Operim*1000+MaskA*50));
    imwrite(ShowOuter,[TPN 'images\ShowOuter.tif'],'Compression','none')
    
    Near=(DistMap>10) & (DistMap<(Outer/2));
    Far=(DistMap>(Outer/2));
    
    Idot=sum(sum((DotsA(Near))))/sum(sum(Area(Near)));
    Odot=sum(sum((DotsA(Far))))/sum(sum(Area(Far)));
    Idend=sum(sum((DendA(Near))))/sum(sum(Area(Near)));
    Odend=sum(sum((DendA(Far))))/sum(sum(Area(Far)));
    
    Grad.Idot=Idot;
    Grad.Odot=Odot;
    Grad.Idend=Idend;
    Grad.Odend=Odend;
    
    Grad.IDD=Idot/Idend;
    Grad.ODD=Odot/Odend;
    Grad.RDD=Grad.IDD/Grad.ODD;
    
    %save([TPN 'Grad.mat'],'Grad')
    
end
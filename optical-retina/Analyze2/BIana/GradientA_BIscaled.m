%%Get Numbers for Inner Half, Outer Half, Inner Third, Outer Third, 
%%inner half - CB, Outer Half / Inner Half - CB

clear all
'Get Cells'
load(['UseCells.mat'])
load('Res.mat')
[Age ONs OFFs BIs Others]=GetCellInfo(UseCells);
GetBIs=find(BIs>0);
BIAge=Age(BIs>0);
'Cells Gotten'
%}


c=0;

NumCells=size(GetBIs,1);

II=.1; %inner boarder

for k = 1:NumCells
    k
    clear  Area Dots Dend DistMap Outer Mask Sum In Out Grad
    TPN = char(UseCells(GetBIs(k),:)); 
    
    
    load([TPN 'UseSc.mat'])
    Cent=Use.Cent;
    clear Use
    
    %% Run from areas
    load([TPN 'CAbiaSC.mat'])
    clear Area Dots Dend
        Area(:,:,1)=CA.Arbor(1).Territory ;
        Area(:,:,2)=CA.Arbor(2).Territory ;
        Dots(:,:,1)=CA.Arbor(1).DotMap;
        Dots(:,:,2)=CA.Arbor(2).DotMap;
        Dend(:,:,1)=CA.Arbor(1).DendMap ;
        Dend(:,:,2)=CA.Arbor(2).DendMap ;
        %image(Dend(:,:,1)*10),pause
           
    

%% Find arbor outer
    DistMap=zeros(size(Area,1),size(Area,2));
    for y = 1:size(Area,1)
        for x = 1: size(Area,2)
            DistMap(y,x)=sqrt((y-Cent(1))^2+(x-Cent(2))^2);
        end
    end
    
clear Outer
for a = 1:2 %run arbors seperately
    TotDend=sum(sum(Dend(:,:,a)));    
    for d = 10:max(DistMap(:))
        Mask=DistMap(:,:)<=d;
        Sum=sum(sum(Dend(:,:,a).*Mask));
        if Sum >=(TotDend * 0.98)
            Outer(a)=d
            break
        end %end if rad big enough
    
    end %end run all dists
    
end %end run both arbors


    
    
%% Chart all
    ilim=10; %inner limit of search
    for a = 1:2
        In(:,:,a)=DistMap>10 & DistMap<(Outer(a)/2);
        Out(:,:,a)=DistMap>(Outer(a)/2);
    end    
    OArea=shiftdim(sum(sum(Area.*Out,1),2),2);
    IArea=shiftdim(sum(sum(Area.*In,1),2));
    Idot=shiftdim(sum(sum(Dots.*In,1),2));
    Odot=shiftdim(sum(sum(Dots.*Out,1),2));
    Idend=shiftdim(sum(sum(Dend.*In,1),2));
    Odend=shiftdim(sum(sum(Dend.*Out,1),2));
    IDD=Idot./Idend;
    ODD=Odot./Odend;
    IdotA=Idot./IArea;
    OdotA=Odot./OArea;
    IdendA=Idend./IArea;
    OdendA=Odend./OArea;
    
    
    dDD=(IDD-ODD)./((IDD+ODD)); %~~~~!!!!!!!!!! Scale?
    dDot=(Idot-Odot)./((Idot+Odot));
    dDend=(Idend-Odend)./((Idend+Odend));
    
    Grad.Idot=Idot;
    Grad.Odot=Odot;
    Grad.Idend=Idend;
    Grad.Odend=Odend;
    Grad.dDD=dDD;
    Grad.dDot=dDot;
    Grad.dDend=dDend;
    Grad.Outer=Outer
    
    Grad.IDD=IDD;
    Grad.ODD=ODD;
    Grad.RDD=IDD./ODD;
    Grad.Rdots=IdotA./OdotA;
    Grad.Rdend=IdendA./OdendA;
    
    
    Grad
    save([TPN 'GradABiSc.mat'],'Grad')
    

%%
    
    
    
    
end
clear all
load('UseCells.mat')
load('Res.mat')
[Age ONs OFFs BIs Others]=GetCellInfo(UseCells);
FindBi=find(BIs)

c=0;
clear Area AreaSc
for k = 1:size(FindBi,1)
    TPN = char(UseCells(FindBi(k)));
    c=c+1;


    %    load([TPN 'CA.mat'])
    %    Area(k,1)=CA.Area;
    %    load([TPN 'Use.mat'])
    %    TotL(k,1)=sum(Use.Length);
    %    load([TPN 'CASc.mat'])
    %    AreaSc(k,1)=CA.Area;
    %    load([TPN 'UseSc.mat'])
    %    TotLSc(k,1)=sum(Use.Length);

    load([TPN 'CABia.mat'])
    Area(k,:)=CA.Area;
    TotDend(c,1)=sum(sum(CA.Arbor(1).DendMap));
    TotDend(c,2)=sum(sum(CA.Arbor(2).DendMap));
    TotDot(c,1)=sum(sum(CA.Arbor(1).DotMap));
    TotDot(c,2)=sum(sum(CA.Arbor(2).DotMap));
    DendMap=CA.Arbor(1).DendMap;
    
    load([TPN 'CAbiaSc.mat'])
    AreaSc(k,:)=CA.Area;
    TotDendSc(c,1)=sum(sum(CA.Arbor(1).DendMap));
    TotDendSc(c,2)=sum(sum(CA.Arbor(2).DendMap));
    TotDotSc(c,1)=sum(sum(CA.Arbor(1).DotMap));
    TotDotSc(c,2)=sum(sum(CA.Arbor(2).DotMap));
    DendMapSc=CA.Arbor(1).DendMap;    
  
    load([TPN 'ONOFFa.mat'])
    load([TPN 'Use.mat'])
    TotL(k,:)=sum(Use.Length);
   
    DotNumA(c,1)=sum(ONOFF.DotLayer==1);
    DotNumA(c,2)=sum(ONOFF.DotLayer==2);
    TotLengthA(c,1)=sum(Use.Length(ONOFF.MidLayer==1));
    TotLengthA(c,2)=sum(Use.Length(ONOFF.MidLayer==2));
        
    load([TPN 'UseSc.mat'])
    TotLSc(k,1)=sum(Use.Length);
    Resolution(c)=Use.Res;

    if Resolution(c) < 0.104
        image((DendMapSc-DendMap)*100 +100)
        %pause
    end
    
    
    DotNumASc(c,1)=sum(ONOFF.DotLayer==1);
    DotNumASc(c,2)=sum(ONOFF.DotLayer==2);
    TotLengthASc(c,1)=sum(Use.Length(ONOFF.MidLayer==1));
    TotLengthASc(c,2)=sum(Use.Length(ONOFF.MidLayer==2));


    load([TPN 'GradABiSc'])
    RDDSc(c,:)=Grad.RDD;
    ODDSc(c,:)=Grad.ODD;
    RdotSc(c,:)=Grad.Rdots;
    RdendSc(c,:)=Grad.Rdend;
    IdotSc(c,:)=Grad.Idot;
    OdotSc(c,:)=Grad.Odot;
    IdendSc(c,:)=Grad.Idend;
    OdendSc(c,:)=Grad.Odend;
    OuterSc(c,:)=Grad.Outer;

    load([TPN 'GradABi'])
    RDD(c,:)=Grad.RDD;
    ODD(c,:)=Grad.ODD;
    Rdot(c,:)=Grad.Rdots;
    Rdend(c,:)=Grad.Rdend;
    Idot(c,:)=Grad.Idot;
    Odot(c,:)=Grad.Odot;
    Idend(c,:)=Grad.Idend;
    Odend(c,:)=Grad.Odend;
    Outer(c,:)=Grad.Outer;

    k
end

scatter(Area(:,1),AreaSc(:,1))
scatter(Area(:,2),AreaSc(:,2))
scatter(Idend(:,1),IdendSc(:,1))
scatter(Odend(:,1),OdendSc(:,1))
scatter(Idot(:,1),IdotSc(:,1))
ylim([0 2500])
xlim([0 2500])

scatter(Outer,OuterSc)


scatter(TotDend(:,1),TotDendSc(:,1))
scatter(TotLengthA(:,2),TotLengthASc(:,2))

scatter(TotDend(:,1),TotLengthA(:,1))
scatter(TotDendSc(:,1),TotLengthASc(:,1))

Area2=Area.* (Res(:,2).^2/0.103^2)

%Percent Difference
PercentDifference=(AreaSc-Area2)./Area2
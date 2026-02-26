clear all
load('UseCells.mat')
load('Res.mat')


c=0;
clear Area AreaSc
for k = 1:size(UseCells,1)
   TPN = char(UseCells(k));
   c=c+1;
   
   
%    load([TPN 'CA.mat'])
%    Area(k,1)=CA.Area;
%    load([TPN 'Use.mat'])
%    TotL(k,1)=sum(Use.Length);
%    load([TPN 'CASc.mat'])
%    AreaSc(k,1)=CA.Area; 
%    load([TPN 'UseSc.mat'])
%    TotLSc(k,1)=sum(Use.Length);

   load([TPN 'CA.mat'])
   Area(k,:)=CA.Area;
   %TotDend=sum(sum(CA.Arbor(1).DendMap))
   load([TPN 'Use.mat'])
   
   TotL(k,:)=sum(Use.Length);
   load([TPN 'CASc.mat'])
   AreaSc(k,:)=CA.Area; 
   load([TPN 'UseSc.mat'])
   TotLSc(k,1)=sum(Use.Length);


       load([TPN 'UseSc.mat'])
       Res(c)=Use.Res;
   
       load([TPN 'GradASc'])
       RDDSc(c,:)=Grad.RDD;
       ODDSc(c,:)=Grad.ODD;

       IdotSc(c,:)=Grad.Idot;
       OdotSc(c,:)=Grad.Odot;
       IdendSc(c,:)=Grad.Idend;
       OdendSc(c,:)=Grad.Odend;
       OuterSc(c,:)=Grad.Outer; 
       
       
       load([TPN 'GradA'])
       RDD(c,:)=Grad.RDD;
       ODD(c,:)=Grad.ODD;
       Idot(c,:)=Grad.Idot;
       Odot(c,:)=Grad.Odot;
       Idend(c,:)=Grad.Idend;
       Odend(c,:)=Grad.Odend;
       Outer(c,:)=Grad.Outer;
   
   
   k
end

scatter(Idend./Odend,IdendSc./OdendSc)
scatter(Outer,OuterSc)
Area2=Area.* (Res(:,2).^2/0.103^2)

%Percent Difference
PercentDifference=(AreaSc-Area2)./Area2
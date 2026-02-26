%%simulate sampling segments with different dot densities. 

clear all

SkelLength=300;
AveDens=.1;
Interval=1/AveDens;
NumDot=fix(SkelLength*AveDens);
RandSkel=rand(NumDot,1)*SkelLength;
EvenSkel=(Interval/2:Interval:SkelLength-Interval/2)';
ClustSkel=EvenSkel(2:3:size(EvenSkel,1));
ClustSkel=[ClustSkel ;ClustSkel ;ClustSkel ];
Scat=randn(size(ClustSkel,1),1)*(Interval/2); %list of random with normal distribution multiplied times standard dev
ClustSkel=ClustSkel+Scat;


UseSkel=RandSkel;

%{
ShowSkel=zeros(SkelLength,1);
for i = 1 : NumDot
    ShowSkel(fix(UseSkel(i))+1)=ShowSkel(fix(UseSkel(i))+1)+1;
end
image(ShowSkel*50)
%}

Samp=1000;
clear Hits
SegLength=Interval;
for i = 1:Samp
   Seg=rand*(SkelLength-Interval);
   Hits(i)=sum((UseSkel>=Seg)&(UseSkel<=(Seg+Interval)));
    
end
hist(Hits,0:1:6)

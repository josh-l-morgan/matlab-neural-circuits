function [DepthDev] = axStrat(TPN,xyum,zum)
%
% clear all
%
% TPN = GetMyDir;
% yxum=0.103;
% zum=0.3;
%


LookRad=30;

load([TPN 'data\AllSeg.mat'])
Mids=mean(AllSeg,3);

load([TPN 'DotSelect.mat'])

DotPos=DotSelect.Pos;
DPos(:,1:2)=DotPos(:,1:2)*xyum;
DPos(:,3)=DotPos(:,3)*zum;


clear DepthDev

NumDots=zeros(size(DPos,1),1);
SDevDots=NumDots;

for i = 1: size(DPos,1)
    Dist=sqrt((DPos(:,1)-DPos(i,1)).^2 + (DPos(:,2)-DPos(i,2)).^2);
    Near=Dist<LookRad;
    Nears=DPos(Near,:);
    NumDots(i)=size(Nears,1);
    SDevsDots(i)=std(Nears(:,3));
end

'Find mids over area'
NumMids=zeros(size(AllSeg,1),1);
SDevsMids=NumMids;
for i = 1: size(y,1)
    Dist=sqrt((Mids(:,1)-y(i)).^2 + (Mids(:,2)-x(i)).^2);
    Near=Dist<LookRad;
    Nears=Mids(Near,:);
    NumMids(i)=size(Nears,1);
    SDevsMids(i)=std(Nears(:,3));
end

OKs=(NumDots>10) & (NumMids > 40);
DepthDev.DotsA=mean(SDevsDots(OKs));
DepthDev.DendA=mean(SDevsMids(OKs));
DepthDev;



'Find Each Dot'
%Run from Dot perspective
NumDots=zeros(size(DPos,1),1);
DtoDots=NumDots;
for i = 1:size(DPos,1)
    Dist=sqrt((DPos(:,1)-DPos(i,1)).^2 + (DPos(:,2)-DPos(i,2)).^2);
    Near=Dist<LookRad;
    Nears=DPos(Near,:);
    NumDots(i)=size(Nears,1);
    meanDepth=median(Nears(:,3));
    DtoDots(i)=abs(meanDepth-DPos(i,3));
end

'Find Each Mid'
%Run from Dot perspective
NumMids=zeros(size(Mids,1),1);
DtoMids=NumMids;
for i = 1:size(Mids,1)
    Dist=sqrt((Mids(:,1)-Mids(i,1)).^2 + (Mids(:,2)-Mids(i,2)).^2);
    Near=Dist<LookRad;
    Nears=Mids(Near,:);
    NumMids(i)=size(Nears,1);
    meanDepth=median(Nears(:,3));
    DtoMids(i)=abs(meanDepth-Mids(i,3));
end




DepthDev.DtoDots=mean(DtoDots(NumDots>10));
DepthDev.DtoMids=mean(DtoMids(NumMids>40));
save([TPN 'data\DepthDevNoCBmedian.mat'],'DepthDev')
DepthDev


%% Draw Data
%
%    for i = 1: size(Mids,1)
%       DrawMids(fix(Mids(i,1))+1,fix(Mids(i,2))+1,fix(Mids(i,3))+1)=1;
%    end
%    imwriteNp(TPN,DrawMids,'DrawMids')
%    %}


%% Compare ON vs OFF
%%different from previous ON OFF in that the only cut off is the dividing
%%line between arbors.

function [ONOFF] = CompairONOFF_2008(TPN, Use)


clear ONOFF

%Get Rotation if it exists
if exist([TPN 'Rot.mat']),
    load([TPN 'Rot.mat'])
    DPos=CoRot(Use.DPos,Rot); %#ok<NODEF>
    Mids=CoRot(Use.Mids,Rot);
    Cent=CoRot(Use.Cent,Rot);
else
    DPos = Use.DPos;
    Mids = Use.Mids;
    Cent = Use.Cent;
end




%%Find cell position

load([TPN 'data\Results.mat'])
Top=[]; Bottom=[];
for a = 1:size(Results.Arbor,2)
    Top(a)=Results.Arbor(a).Top;
    Bottom(a)=Results.Arbor(a).Bottom;
end


CellZ=Cent(3);
if abs(CellZ-mean(Top))>abs(CellZ-mean(Bottom));
    Flip=1; else Flip=0;
end
% Fr(k)=Flip;

Middle=(Bottom(2)+Top(1))/2;


Layer=zeros(size(Mids,1),1);
Layer((Mids(:,3)<Middle(1)))=1;
if size(Top,2)>1
    Layer((Mids(:,3)>Middle) )=2;
end

DotLayer=zeros(size(DPos,1),1);
DotLayer((DPos(:,3)>=Top(1)) & (DPos(:,3)<Middle(1)))=1;
DotLayer((DPos(:,3)>Middle) )=2;



ONOFF.DotLayer=DotLayer;
ONOFF.MidLayer=Layer;

ONOFF.length(1)=sum(Use.Length(ONOFF.MidLayer==1));
ONOFF.length(2)=sum(Use.Length(ONOFF.MidLayer==2));
ONOFF.dots(1)=sum(ONOFF.DotLayer==1);
ONOFF.dots(2)=sum(ONOFF.DotLayer==2);
ONOFF.flip = Flip;

save([TPN 'ONOFFa.mat'],'ONOFF')

end













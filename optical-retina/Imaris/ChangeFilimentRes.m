TPN = GetMyDir;
load([TPN 'Skel.mat'])

%%FilStats, SegStats

Res=[0.062; 0.103; 0.115];


ImageInfo= inputdlg({'SkelRes','RealRes'},'FixSkelRes',1,{'0.103','0.115'});

OldRes=str2num(ImageInfo{1});
NewRes=str2num(ImageInfo{2});

aXYZ=FilStats.aXYZ;
aEdges=FilStats.aEdges+1;

aXYZ(:,1:2)=aXYZ(:,1:2)*(NewRes/OldRes);


Lengths=sqrt((aXYZ(aEdges(:,1),1)-aXYZ(aEdges(:,2),1)).^2 ...
           + (aXYZ(aEdges(:,1),2)-aXYZ(aEdges(:,2),2)).^2 ...
           + (aXYZ(aEdges(:,1),3)-aXYZ(aEdges(:,2),3)).^2);
       
Length=sum(Lengths)
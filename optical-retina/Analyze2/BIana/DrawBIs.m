%% Compare ON vs OFF

'Get Cells'
load(['UseCells.mat']);
%[Age ONs OFFs BIs Others]=GetCellInfo(UseCells);
GetBIs=find(BIs);
BIAge=Age(BIs);
'Cells Gotten'


for k = 1:size(UseCells,1)
   TPN=char(UseCells(k));
   load([TPN 'Use.mat']);
   aNum(k,1)=size(Use.Top,2);
    
end

% Find Ambiguous cells
Amb=xor(aNum>1,BIs>0);

gAmb=find(Amb);

TwoArbors=find(aNum>1);
%% Draw Mids

for k = 1:size(GetBIs,1)
   clear CellZ Top Bottom Layer Mids MaxMids DrawMids Side
   TPN=char(UseCells(GetBIs(k)))
   load([TPN 'Use.mat'])
   Mids=Use.Mids;
   
   if exist([TPN 'Rot.mat'])
   load([TPN 'Rot.mat'])
   Mids=CoRot(Mids,Rot);
   end
  
    CellZ=Use.Cent(3);
    
    load([TPN 'data\ResultsComp.mat'])
    Top=[]; Bottom=[];
    for a = 1:size(Results.Arbor,2)
        Top(a)=Results.Arbor(a).Top;
        Bottom(a)=Results.Arbor(a).Bottom;
    end
   
    Layer=zeros(size(Mids,1),1)+50;
    
    Layer((Mids(:,3)>=Top(1)) & (Mids(:,3)<=Bottom(1)))=150;
    if size(Top,2)>1
        Layer((Mids(:,3)>Top(2)) & (Mids(:,3)<=Bottom(2)))=250;
    end
   
    
   %%Transform Mids for Drawing  
   minMids=min(Mids);
   for m=1:3; Mids(:,m)=Mids(:,m)-minMids(m); end
    
   Mids=fix(Mids)+1;
   MaxMids=max(Mids);
     
    
    
   DrawMids=zeros([MaxMids],'uint8');
   for i = 1:size(Mids,1)
       DrawMids(Mids(i,1),Mids(i,2),Mids(i,3))=Layer(i);
   end
   imwriteNp(TPN,DrawMids,'3DMidsR')   
   
   Side=max(DrawMids,[],2);
   Side=shiftdim(Side,2);
   image(Side),pause(.3)
    
end

AmbNames=UseCells(gAmb);







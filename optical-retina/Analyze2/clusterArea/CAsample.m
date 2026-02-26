
clear all

KPN=GetMyDir
Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;

for k = 1:size(Kdir,1)
    TPN = [KPN '\' Kdir(k).name '\']; 
     clear CA BiDotDist DotShow DotFilt DendFilt DotMap DendMap DPos 
     clear Cell Red Green Blue TerFill TerFilt TerPerim TerRay Territory
     clear Mids sMids NN sNN
    if exist([TPN 'find\SG.mat'])
      
    
%% load info
'Load Data'
load([TPN 'Dots.mat'])
load([TPN 'find\SG.mat'])
load([TPN 'data\AllSegCut.mat'])
load([TPN 'data\Results.mat'])
load('.\cmap.mat')

set(gcf,'Colormap',cmap)

%%Extract Dot Positions
OK=SG.passF;
DPos=Dots.Pos(OK,:);
DPos(:,1:2)=DPos(:,1:2)*.103; DPos(:,3)=DPos(:,3)*.3;
Cent=Dots.Im.CBpos;
Cent=Cent*.103;
clear Dots SG


%%Extract Dend positions 
Mids=mean(AllSegCut,3);
Length=sqrt((AllSegCut(:,1,1)-AllSegCut(:,1,2)).^2 ...
           + (AllSegCut(:,2,1)-AllSegCut(:,2,2)).^2 ...
           + (AllSegCut(:,3,1)-AllSegCut(:,3,2)).^2);
clear AllSegCut

%%Extract Depth restrictions
clear Top Bottom
for i = 1: size(Results.Arbor,2)
    Top(i)=Results.Arbor(i).Top;
    Bottom(i)=Results.Arbor(i).Bottom;
end
clear Results

%% create Nearest Node list (nearest node for each dot and distance)
for i = 1:size(DPos,1)
            Ndist=dist(Mids,DPos(i,:)); %find dist from dot to all nodes
            Near=min(Ndist); %find shortest distance
            Nearest=find(Ndist==Near,1); %get node at that distance
            NN(i,:)=Mids(Nearest,:); %assign that node to NearestNode list for dots
            DotToNN(i,:)=Near; %record that distance for posterity  
end

% Get image properties
    minDPos=min(NN); minMids=min(Mids);
    minBoth=fix(min(minDPos,minMids))-1;
    
    %%insure no negs

    for d = 1:2, 
        if minBoth(d)<1;
            NN(:,d)=NN(:,d)-minBoth(:,d); 
            Mids(:,d)=Mids(:,d)-minBoth(:,d); 
        end
    end
    Msize=fix(max(max(NN),max(Mids)))+2;

%%
for a = 1:size(Top,2)+1  % Run Arbors
    
    clear sMids sNN
    if a>1
        sMids=Mids(Mids(:,3)>Top(a-1) & Mids(:,3)<Bottom(a-1),:);
        sNN=NN(NN(:,3)>Top(a-1) & NN(:,3)<Bottom(a-1),:);
        
    else
        sMids=Mids;
        sNN=NN;
    end
        
    sMids=round(sMids);
    sNN=round(sNN);

    %% Draw maps
    'Draw Maps'

    %%Draw Dot map
    DotMap=zeros(Msize(1),Msize(2));
    for i=1:size(sNN,1)
        DotMap(sNN(i,1),sNN(i,2))=DotMap(sNN(i,1),sNN(i,2))+1;
    end
    
    %%Draw Dend map
    DendMap=zeros(Msize(1),Msize(2));
    for i=1:size(sMids,1)
        DendMap(sMids(i,1),sMids(i,2))=DendMap(sMids(i,1),sMids(i,2))+Length(i);
    end

%% Find Territory
    
    %% make distance filter

    AreaS=10;
    Disk=fspecial('disk',AreaS);
    DotFilt=imfilter(DotMap,Disk,'same');
    DendFilt=imfilter(DendMap,Disk,'same');
    BlankFilt=imfilter(double(DendMap>-1000),Disk,'same');
    DotFilt=DotFilt./BlankFilt;
    DendFilt=DendFilt./BlankFilt;

    DendSkel=DendFilt;
    DendSkel(DendMap>0)=0;
    image(DendSkel*1000)

    %DDSkel=DotFilt./DendFilt;
    %DDSkel(DendMap>0)=0;
    %image(DDSkel*500)

    DotSkel=DotFilt;
    DotSkel(DendMap>0)=0;
    image(DotSkel*2000)


    %% Find Ray territory
    image(DendMap*100), pause(.3)

    Disk2=fspecial('disk',5);%min(5,min(AreaS/2,1)));
    TerFilt=imfilter(DendMap,Disk2,'same');

    DFlab=bwlabel(TerFilt);
    for i = 1:max(DFlab(:))
        lSize(i)=size(find(DFlab==i),1); %#ok<AGROW>
    end
    TerFilt=DFlab*0;
    TerFilt(DFlab==find(lSize==max(lSize)))=1;
    image(TerFilt*200),pause(.3)
    clear DFlab lSize
    
    %Close
    Csize=round(AreaS/2);
    SE=strel('disk',Csize);
    Buf=Csize*2;
    [tys txs] = size(TerFilt);
    BufT=zeros(tys+2*Buf,txs+2*Buf);
    BufT(Buf+1:Buf+tys,Buf+1:Buf+txs)=TerFilt;
    BufC=imclose(BufT,SE);
    TerC=BufC(Buf+1:Buf+tys,Buf+1:Buf+txs);
    image((TerFilt+TerC)*100), pause(.3)
    clear Buf BufT BufC
    
    %Remove holes
    TerFill=TerC;
    Frame=TerC*0; Frame(1:tys,1)=1; Frame(1:tys,txs)=1; Frame(1,1:txs)=1; Frame(tys,1:txs)=1;
    TerHole=bwlabel(~TerC);
    for i=1:max(TerHole(:))
        Surround=sum(sum((TerHole==i) .* Frame));
        if ~Surround
            TerFill(TerHole==i)=1;
        end        
    end
    image(TerFill*100),pause(.01)
    clear Frame TerHole Surround 
    
    %Get Perimeter
    TerPerim=bwperim(TerFill,8);
    image((TerFill+TerPerim)*100),pause(.3)

    %DrawRays
    [tpy tpx]=find(TerPerim);
    TerRay=TerFill*0;

    for i = 1: size(tpy,1)
        ydif=Cent(1)-tpy(i); xdif=Cent(2)-tpx(i);
        Long=sqrt(ydif^2 + xdif^2);
        ystep=ydif/Long;
        xstep=xdif/Long;

        for l = 1:.5:Long
            TerRay(round(tpy(i)+ystep*l),round(tpx(i)+xstep*l))=1;
        end
    end
    %image(TerRay*100),pause(.01)
    
    Territory=TerFill;
    DotDist=DotSkel.*Territory;
    image(DotDist*2000)
    
%%   Collect Data

 
    CA.Arbor(a).DotMap=DotMap;
    CA.Arbor(a).DendMap=DendMap;
    CA.Arbor(a).Territory=Territory;
    CA.Arbor(a).DotDist=DotFilt.*Territory;
    CA.Arbor(a).DendDist=DendFilt.*Territory;
    CA.Arbor(a).TerPerim=TerPerim;
    CA.Arbor(a).MeanDotOArea=mean(DotFilt(Territory>0))

end % Run Arbors

%% Manage Information
if size(CA.Arbor,2)>2
    CA.mDoA=mean([CA.Arbor(2).DotDist(CA.Arbor(2).Territory>0) ;...
        CA.Arbor(3).DotDist(CA.Arbor(3).Territory>0)]);
    CA.Area=sum(CA.Arbor(2).Territory(:))+sum(CA.Arbor(3).Territory(:));
    ScaleLamina=CA.Arbor(2).Territory+CA.Arbor(3).Territory;
    ScaleLamina(ScaleLamina==2)=.5;
    CA.Im.DotDist=CA.Arbor(1).DotDist.*ScaleLamina;
    CA.Im.DendDist=CA.Arbor(1).DendDist.*ScaleLamina;
   
else
    CA.mDoA=CA.Arbor(1).MeanDotOArea;
    CA.Area=sum(CA.Arbor(1).Territory(:));
    CA.Im.DotDist=CA.Arbor(1).DotDist;
    CA.Im.DendDist=CA.Arbor(1).DendDist;
end



ShowCell=(CA.Im.DotDist).* ~CA.Arbor(1).DendMap;
ShowCell=ShowCell*1000+TerPerim;
ShowCell=uint8(ShowCell);
image(ShowCell),pause(.3)
CA.Im.ShowCell=ShowCell;

%%Write images
clear Cell
if exist([TPN 'Cell.mat'])
    load([TPN 'Cell.mat'])
    cName=['P' Cell.Age '_' Cell.Type '_' Kdir(k).name '.tif']
    Name=[KPN 'Images\Area\' cName];
else
    Name=[KPN 'Images\Area\' Kdir(k).name '.tif'];
end
if isempty(find(Name=='?'));
    imwrite(ShowCell,cmap,Name,'Compression','none')
end

if size(CA.Arbor,2)>2
   ScaleRG=max(max(CA.Arbor(3).DotDist(:)),max(CA.Arbor(2).DotDist(:)));
   ScaleRG=255/ScaleRG;
   Red=CA.Arbor(2).DotDist*ScaleRG;
   Red=Red+CA.Arbor(2).TerPerim*75;
   Green=CA.Arbor(3).DotDist*ScaleRG;
   Green=Green+CA.Arbor(3).TerPerim*75;
   Blue=~(CA.Arbor(3).Territory & CA.Arbor(2).Territory);
   BiDotDist(:,:,1)=uint8(Red);
   BiDotDist(:,:,2)=uint8(Green);
   BiDotDist(:,:,3)=uint8(Blue)*0;
   CA.Im.BiDotDist=BiDotDist;
   image(BiDotDist),pause(.3)
   
   if exist('Cell')
    cName=['P' Cell.Age '_' Cell.Type '_' Kdir(k).name '.tif']
    Name=[KPN 'Images\Bi\' cName];
    else
        Name=[KPN 'Images\Bi\' Kdir(k).name '.tif'];
    end
    if isempty(find(Name=='?'));
        imwrite(BiDotDist,Name,'Compression','none')
    end
end
    AreaResults(k)=CA.mDoA;
    Age(k)=str2num(Cell.Age);
    save([TPN 'CA.mat'],'CA')
    
    
      [ys xs]=size(DotMap);
      
      DotDist=CA.Im.DotDist;
      DendDist=CA.Im.DendDist;
      DotMap=CA.Arbor(1).DotMap;
      DendMap=CA.Arbor(1).DendMap;
      CA4=zeros(ys*2,xs*2);
      CA4(1:ys,1:xs)=DotDist*2000;
      CA4(1:ys,xs+1:xs*2)=DendDist*600;
      CA4(ys+1:ys*2,1:xs)=((DotMap>0)+(DendMap>0))*100;
      CA4(ys+1:ys*2,xs+1:xs*2)=DotDist./DendDist*400;
      image(CA4),pause(.3)
      
      if exist([TPN 'Cell.mat'])
            load([TPN 'Cell.mat'])
            cName=['P' Cell.Age '_' Cell.Type '_' Kdir(k).name '.tif']
            Name=[KPN 'Images\C4\' cName];
        else
            Name=[KPN 'Images\C4\' Kdir(k).name '.tif'];
        end
        if isempty(find(Name=='?'));
            imwrite(CA4,cmap,Name,'Compression','none')
        end
      
    
    

end %if find exists, run cell
pause(.3)
end % Run all cells

scatter(Age,AreaResults)


     
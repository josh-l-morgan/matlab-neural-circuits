%{
Published in : Developmental patterning of glutamatergic synapses onto retinal ganglion cells
Josh L Morgan, Timm Schubert and Rachel OL Wong*
Skeletonizes binary mask of dendrites by propogating activated wave
through the mask.  Output is TPN/data/AllSeg.mat and can be replaced by
similar file generated in Imaris.
%}
clear all
subplot(1,1,1)
TPNm = GetMyDir;
dTPNm=dir(TPNm); dTPNm=dTPNm(3:size(dTPNm,1));

slash = regexp(TPNm,'\')
TPN = TPNm(1:slash(end-1));
if ~exist([TPN 'data']), mkdir([TPN 'data']),end
if ~exist([TPN 'temp']), mkdir([TPN 'temp']),end

xyum=4;
zum=30;

pause(.1)


aspect=zum/xyum;% ratio of z to xy dimentions
minObSize= 50; %% minimum size of object to be measured (default 20)
minFillSize = 10; %% Minimum size of continuous wave object (default 4)
maxSegLength = 5;  %% Maximum length of segments (default 2)

%% READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'reading image'


if size(dTPNm,1)==1
    Ir=tiffread2([TPNm dTPNm(1).name]);
    D=zeros(size(Ir(1).data,1),size(Ir(1).data,2),size(Ir,2));
    for i = 1: size(Ir,2)
        D(:,:,i)=Ir(i).data;
    end
    clear Ir
else
    
    clear I Ic
    D(:,:)=imread([TPNm dTPNm(1).name]); %read
    D(1,1,size(dTPNm,1))=0;
    c=0;
    for i=1:size(dTPNm,1)
        nam=dTPNm(i).name;
        naml=length(nam);
        if nam(naml-3:naml)=='.tif'
            c=c+1;
            D(:,:,c)=imread([TPNm nam]);
        end
        PercentRead=i/size(dTPNm,1)*100
    end
    
   D=D==1;
 
end

%Find New sizes
save([TPN 'D.mat'],'D')
[ys,xs,zs]=size(D);


%% Run Skeletonize

'Filling to find Skeleton'

SG=0; %start segment counter
PixDone=0; %Counter for Fill Progress
AllPix=sum(D(:)); %Count all pixels to do
rep=0; %start another damn counter
D=(smooth3(D,'gaussian',5)>=.3); %Thats the images !!!!!!!!!!!!!!!!!!!!>=1?

%%Buffer D
Dbuf=zeros(ys+2,xs+2,zs+2,'uint8');
Dbuf(2:ys+1,2:xs+1,2:zs+1)=D;
D=Dbuf;
clear Dbuf
[ys,xs,zs]=size(D);

sumD=sum(D,3);
colormap gray(256)
image(sumD*15/(mean(sumD(:))))

[D,ob]=bwlabeln(D,6); %Label Mask Objects -Six Connected
tic
for o=1:ob  %Run for All objects
    if sum(sum(sum(D==o)))>minObSize;
    %% Display Current object
    %  image((max(D,[],3)>0)*75 + max(D==o,[],3)*300),pause(.01)
       
    %%Find minimum length as princomp
    [vy,vx,vz]=ind2sub(size(D),find(D==o)); %Vectorize object 
   
    clear V
    V(:,1)=vy;V(:,2)=vx;V(:,3)=vz; %make single vector
    clear vy vx vz
    [Co,Sc]=princomp(V);  %find principal components of object
    minV=find(Sc(:,1)==min(Sc(:,1)),1); %find position of min extreme position
    maxV=find(Sc(:,1)==max(Sc(:,1)),1); %find position of max extreme position
    startV=find(abs(Sc(:,1))==max(abs(Sc(:,1))),1); %find most extreme point for later starting point
    minLength=Sc(maxV,1)-Sc(minV,1); %Find minimum Length in Voxels
    clear Co Sc
    
    %%Start Filling
    con6=[-1 0 0; 1 0 0; 0 -1 0 ; 0 1 0; 0 0 -1; 0 0 1]; %6con matrix
    NewF=zeros(1000,3);
    NewF(1,:)=V(startV,:); %make starting point
    Sold=V(startV,:); %make Skeleton Starting Point
        
    Fnew=D*0;
    count=1;
   while count>0  %while there are new voxels to do
        
    PrevF=NewF(NewF(:,1)>0,:); NewF=NewF*0;  %update list of new fill points
    Fnew=Fnew*0; %Create matrix for new pixels
    count=0;  %reset counter
    for f=1:size(PrevF,1) %run all Previously filled points
    
    %%transform    
    Trans=[PrevF(f,1)+con6(:,1) PrevF(f,2)+con6(:,2) PrevF(f,3)+con6(:,3)];
    
    for t = 1:size(Trans,1)
       Grab(t)=D(Trans(t,1),Trans(t,2),Trans(t,3));
       D(Trans(t,1),Trans(t,2),Trans(t,3))=0;
       Fnew(Trans(t,1),Trans(t,2),Trans(t,3))=1;
    end
    
    
%     Ids=sub2ind([ys xs zs],Trans(:,1),Trans(:,2),Trans(:,3));
%     Grab=D(Ids)>0;    
%     D(Ids)=0; %Fill spot in registry
%     Fnew(Ids)=1;

    GrabF=find(Grab);
    GrabFL=length(GrabF);
    PixDone=PixDone+GrabFL;
    NewF(count+1:count+GrabFL,1:3)=Trans(Grab>0,:);
    count=count+GrabFL;
    
    end % Run new points
            
    %Check new fill objects
    [Fnew,things]=bwlabeln(Fnew,26); %Label new fill objects (with 26 connectivity)
    for t=1:things %Do all objects
        if sum(sum(sum(Fnew==t)))>minFillSize %minimum centroid size
            [ly,lx,lz]=ind2sub(size(Fnew),find(Fnew==t)); %Find position of all pixels in object
            my=mean(ly); mx=mean(lx); mz=mean(lz); %Find centroid
            Snew(t,1)=my;Snew(t,2)=mx; Snew(t,3)=mz; %record new points
            SDist=sqrt((Sold(:,1)-my).^2+ (Sold(:,2)-mx).^2 + (Sold(:,3)-mz).^2); %find voxel distance from new to all old
            source=find(SDist==min(SDist),1);

            %Write Segments
            SG=SG+1;
            Seg(SG,1,1)=Sold(source,1);Seg(SG,2,1)=Sold(source,2); Seg(SG,3,1)=Sold(source,3);
            Seg(SG,1,2)=my; Seg(SG,2,2)=mx; Seg(SG,3,2)=mz;
            Obj(SG)=o; %record objects
            
        end %end if big enough to be worth new centroid
    end
    if exist('Snew','var'), Sold=Snew; clear Snew; end  %update skeleton points list: new becomes old
   
     if mod(SG,30)==0,  image((sumD>0)*255 - (sum(D,3)>0)*100),pause(.01),end
    
   end % Repeat fill
    

   end %if object is big enough   
end %run all objects

clear AllPix D Fnew

Duration = toc


%% Convert to real distances
if exist('Seg','var') % If any Segments were found
%%Unbuffer
Seg=Seg-1;
clear SU
SU(:,1:2,:)=Seg(:,1:2,:)*xyum;
SU(:,3,:)=Seg(:,3,:)*zum;
clear Seg
save([TPN 'temp/SU.mat'],'SU')

%%Remove overly long segments
SegLength=sqrt((SU(:,1,1)-SU(:,1,2)).^2 + (SU(:,2,1)-SU(:,2,2)).^2 + (SU(:,3,1)-SU(:,3,2)).^2);
SU=SU(SegLength<=maxSegLength,:,:);



'Skeleton Complete'

%% Find Nodes
'Finding Nodes'


%% Map Nodes

%%Find Connections
SkelSize=size(SU,1);
for i=1:SkelSize
   
    %%Check how many times each tip is an origin
    TipFind=abs(SU(:,1,1)-SU(i,1,2))+abs(SU(:,2,1)-SU(i,2,2))+abs(SU(:,3,1)-SU(i,3,2));
    Con(i,2)=sum(TipFind==0); %add occurances as origins
    TipFind=abs(SU(:,1,2)-SU(i,1,2))+abs(SU(:,2,2)-SU(i,2,2))+abs(SU(:,3,2)-SU(i,3,2));
    Con(i,2)=Con(i,2)+sum(TipFind==0); %add occurances as tips
    
    %%Check how many times each origin is an origin
    OrigFind=abs(SU(:,1,1)-SU(i,1,1))+abs(SU(:,2,1)-SU(i,2,1))+abs(SU(:,3,1)-SU(i,3,1));
    Con(i,1)=sum(OrigFind==0);
    OrigFind=abs(SU(:,1,2)-SU(i,1,1))+abs(SU(:,2,2)-SU(i,2,1))+abs(SU(:,3,2)-SU(i,3,1));
    Con(i,1)=Con(i,1)+sum(OrigFind==0);
    Con(i,3)=sum(OrigFind==0);  %Secret Origin (0 if never a tip)
end




%% Remove some Intermediates
'Removing intermediates'
n=0; %node counter
minDist=.5; %assign minimum distance
clear nodes Object
for i=1:SkelSize
   %% Search list for Origin points
   if Obj(i)==29, 'running 29', end
    
    if Con(i,1)>2 | Con(i,3)==0,  Origin=SU(i,:,1); %if current origin is branch point or never tip
        if Obj(i)==29, '29 passed Origin Test', end
        Link=Origin;  %make the origin the first link
        for f=i:SkelSize  %run every subsequent segment from current point to end of list. 
             if sum(abs(SU(f,:,1)-Link))==0 %Look for the next segment
                 if Obj(i)==29, 'Found next Link', end
                Tip=SU(f,:,2); %if the origin is a link, then assign the new tip
                Link=Tip;
                Dist=sqrt((Tip(1)-Origin(1))^2 + (Tip(2)-Origin(2))^2 + (Tip(3)-Origin(3))^2); %find distance
                 %%Condition 1 =  Tip is end point
                if Con(f,2) ~=2 %If tip is end point or branch point
                   n=n+1;
                   nodes(n,:,1)=Origin; nodes(n,:,2)=Tip; %Then assign segment to nodes
                   Object(n)=Obj(i); % record source Object
                   break %leave f to look down list for next origin
                   if Obj(i)==29, 'Hit end point and wrote segment', end
                elseif Dist>minDist %If far enough away to be segment
                   n=n+1;
                   nodes(n,:,1)=Origin; nodes(n,:,2)=Tip; %Then assign segment to nodes
                   Object(n)=Obj(i);  %Record Source Object
                   Origin=Tip; %update tip as new origin
                   if Obj(i)==29, 'Hit min dist and wrote Segment', end
                end  %Tip condition      
             end %if its the right segment
                
         end %if current segment starts with the link

    end %End if i is an origin.
end


clear Con Dist SU


%% ReMap Nodes

%%Find Connections
clear TipFind OrigFind Con2
NodeSize=size(nodes,1);
for i=1:NodeSize
   
    %%Check how many times each tip is an origin
    TipFind=abs(nodes(:,1,1)-nodes(i,1,2))+abs(nodes(:,2,1)-nodes(i,2,2))+abs(nodes(:,3,1)-nodes(i,3,2));
    Con2(i,2)=sum(TipFind==0); %add occurances as origins
    TipFind=abs(nodes(:,1,2)-nodes(i,1,2))+abs(nodes(:,2,2)-nodes(i,2,2))+abs(nodes(:,3,2)-nodes(i,3,2));
    Con2(i,2)=Con2(i,2)+sum(TipFind==0); %add occurances as tips
    
    %%Check how many times each origin is an origin
    OrigFind=abs(nodes(:,1,1)-nodes(i,1,1))+abs(nodes(:,2,1)-nodes(i,2,1))+abs(nodes(:,3,1)-nodes(i,3,1));
    Con2(i,1)=sum(OrigFind==0);
    OrigFind=abs(nodes(:,1,2)-nodes(i,1,1))+abs(nodes(:,2,2)-nodes(i,2,1))+abs(nodes(:,3,2)-nodes(i,3,1));
    Con2(i,1)=Con2(i,1)+sum(OrigFind==0);
    Con2(i,3)=sum(OrigFind==0);  %Secret Origin (0 if never a tip)
end



% %% Connect the Unconnected
% 'Connecting Unconnected'
% maxbridge=3; %maximum distance to bridge a gap in microns
% [Ends Type]=find(Con2(:,1:2)==1); Es=size(Ends); %Find all one connected points
% b=0; %reset bridge counter
% Bridged=0; %matrix for recording linked objects
% clear Dist bridge
% for i=1:size(Ends,1)
%     %%Find distance to all points
%     Dist=sqrt((nodes(:,1,:)-nodes(Ends(i),1,Type(i))).^2+ (nodes(:,2,:)-nodes(Ends(i),2,Type(2))).^2 + (nodes(:,3,:)-nodes(Ends(i),3,Type(i))).^2); 
%     home=Object(Ends(i));  %Identify home object of tip
%     Dist(Object'==home,1,:)=2*maxbridge; %eliminate segs of same object by increasing distance
%     minDist1=min(Dist(:,1,1));
%     minDist2=min(Dist(:,1,2));
%     if minDist1<minDist2 & minDist1<maxbridge
%        Btarget=find(Dist(:,1,1)==minDist1,1); %find nearest ID
%           b=b+1; %increace counter
%           bridge(b,:,1)=nodes(Btarget,:,1);  %make nearest the new origin
%           bridge(b,:,2)=nodes(Ends(i),:,Type(i)); %make current end a new tip
%           Bridged(b,1)=home; Bridged(b,2)=Object(Btarget); %record Briged objects
%      elseif minDist1 <= minDist2 & minDist2<maxbridge
%        Btarget=find(Dist(:,1,2)==minDist2,1);
%            b=b+1;
%            bridge(b,:,1)=nodes(Btarget,:,2); %make nearest the new origin
%            bridge(b,:,2)=nodes(Ends(i),:,Type(i)); %make current end a new tip
%            Bridged(b,1)=home; Bridged(b,2)=Object(Btarget); %record Briged objects
%     end
% end %End i, search all tips
% 
% %%Eliminate extra bridges
% if exist('bridge','var') %if bridge exists
% 
% 
% %%Assign Bridge IDs
% id=0; %start ID counter
% bid=zeros(size(bridge,1),1); %create list for bridge IDs
% for i=1:size(bridge,1) %run all bridges
%     if bid(i)==0 % if starting bridge hasnt been IDed
%         id=id+1; %increace id counter
%     for j=1:size(bridge,1)  %check against all bridges
%         if bid(j)==0  %if target hasnt been IDed
%             dif1=sum(abs(Bridged(j,1)-Bridged(i,1))+abs(Bridged(j,2)-Bridged(i,2))); %check if same connection
%             dif2=sum(abs(Bridged(j,1)-Bridged(i,2))+abs(Bridged(j,2)-Bridged(i,1))); %check if reverse connection
%             if dif1 ==0 | dif2==0, bid(j)=id; end %if its the same connection, assign current ID
%         end %if j IDed
%     end %j, run all targets
%     end %if i has been IDed
% end %run all bridges
% 
% %%find Shortest Bridges
% Bdist=sqrt((bridge(:,1,1)-bridge(:,1,2)).^2+(bridge(:,2,1)-bridge(:,2,2)).^2+(bridge(:,3,1)-bridge(:,3,2)).^2); %find all distances
% bridges=zeros(max(bid),3,2);
% for i=1:max(bid)
%     testdist=Bdist; % make temp bridge list
%     testdist(bid~=i)=maxbridge*2; %eliminate other ids
%     BmTarget=find(testdist==min(testdist),1); %find shortest bridge
%     bridges(i,:,:)=bridge(BmTarget,:,:); %enter shortest bridge
% end
% 
% bridge=bridges; %replace redundant bridge list with pruned bridge list
% end %if bridge exists
% 
% clear Bdist BmTarget Bridged Btarget Dist Ends bridges

%% Create Final vector map
if exist('bridge','var')
    NewSeg=[nodes ; bridge];  %Create variable with all identified segments
else 
    NewSeg=nodes
end


%%Diagnostic, get rid of later
Sc=1/xyum;


end %if any segments were found

%% Cheap
AllSeg=NewSeg;


%%move AllSeg from temp to Data
if AllSeg(1,1,1)==0, AllSeg=AllSeg(2:size(AllSeg,1),:,:); end %Get rid of any spacer
save([TPN 'data/AllSeg.mat'],'AllSeg')
% 
% SegLength=sqrt((AllSeg(:,1,1)-AllSeg(:,1,2)).^2 +(AllSeg(:,2,1)-AllSeg(:,2,2)).^2 +(AllSeg(:,3,1)-AllSeg(:,3,2)).^2);
% save([TPN 'data/SegLength.mat'],'SegLength')
% TotalLength=sum(SegLength);
% save([TPN 'data/TotalLength.mat'],'TotalLength')


%% Finish
% 
% [TPN(size(TPN,2)-6:size(TPN,2)-1)]
% SkeletonizedAt=uint16(clock)
% save([TPN 'data/SkeletonizedAt.mat'],'SkeletonizedAt')

'Done Skeletonizing'

%% Data Drawing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw Nodes
% 'Drawing Data'
% load([TPN 'data/AllSeg.mat'])
% Sc=1/xyum; %scalling factor to decide image size each pix width in um = 1/Sc
% AllSegI=0;
% 
% for i=1:size(AllSeg,1)
%     AllSegI(round(AllSeg(i,1,1)*Sc)+1,round(AllSeg(i,2,1)*Sc)+1,round(AllSeg(i,3,1)*Sc)+1)=1;  %Draw all origins.
%     AllSegI(round(AllSeg(i,1,2)*Sc)+1,round(AllSeg(i,2,2)*Sc)+1,round(AllSeg(i,3,2)*Sc)+1)=1;  %Draw all tips.
% end
% 
% maxAllSegI=(max(AllSegI,[],3)>0);
% image(maxAllSegI*300),pause(2)
% imwriteNp(TPN,AllSegI,'Nodes')
% 


% %% Draw Skeleton
% SkelRes=1;
% Length=zeros(size(AllSeg,1),1);
% Sc=1/xyum;
% Skel=zeros(fix(ys*xyum)+1,fix(xs*xyum)+1,fix(zs*zum)+1,'uint8');
% %clear AllSegI
% for i=1:size(AllSeg,1)
%     Dist=sqrt((AllSeg(i,1,1)-AllSeg(i,1,2))^2 + (AllSeg(i,2,1)-AllSeg(i,2,2))^2 + (AllSeg(i,3,1)-AllSeg(i,3,2))^2); %find distance
%     Length(i)=Dist;
%       devs=max(1,round(Dist/SkelRes)); %Find number of subdivisions
%     for d=1:devs+1
%         sy=AllSeg(i,1,1)+((AllSeg(i,1,2)-AllSeg(i,1,1))/devs)*(d-1);
%         sx=AllSeg(i,2,1)+((AllSeg(i,2,2)-AllSeg(i,2,1))/devs)*(d-1);
%         sz=AllSeg(i,3,1)+((AllSeg(i,3,2)-AllSeg(i,3,1))/devs)*(d-1);
%         Skel(round(sy*Sc)+1,round(sx*Sc)+1,round(sz*Sc)+1)=1; %draw Skel
%     end
% end
% clear Dist
% maxSkel=(max(Skel,[],3)>0);
% % imwriteNp(TPN,Skel,'Skel')

%% colorize skeleton
% Skel=uint8(Skel);
% for i=1:size(Skel,3)
%     Dmap=uint8(Skel(:,:,i)>0);
%     Dmap(Dmap>0)=i;
%     Skel(:,:,i)=Dmap;
% end
% 
% imwriteNp(TPN,Skel,'SkelDepth')


%% Image Isosurface
% 
% isosurface(Skel,0), axis equal, view(3)
% camlight, lighting gouraud, title('Skeleton') 
% clear Bridge 
% 
% Duration
% %clear BlockBuffer BlockSize
% %}

%% NOTES
%{



%}

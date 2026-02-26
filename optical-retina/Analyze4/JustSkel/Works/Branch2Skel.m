load([TPN 'Branches.mat'])

%% Convert Branches to Segs
clear Seg
SG=0;
for b = 1: length(Branch)
    
    GetBranch=Branch{b};
    GetSizes=Sizes{b};
    if sum(GetSizes)>minFillSize
        for n = 2+SkipJoints: size(GetBranch,1)
            SG = SG+1;
            Seg(SG,:,1)=GetBranch(n-1,:);
            Seg(SG,:,2)=GetBranch(n,:);
            Obj(SG)=Source(bn,1);
        end
    end
end


%% Convert to real distances
if exist('Seg','var') % If any Segments were found
    %%Unbuffer
    
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
        NewSeg=nodes;
    end


    %%Diagnostic, get rid of later
    Sc=1/xyum;


end %if any segments were found

%% Cheap
AllSeg=NewSeg;


%%move AllSeg from temp to Data
if AllSeg(1,1,1)==0, AllSeg=AllSeg(2:size(AllSeg,1),:,:); end %Get rid of any spacer
save([TPN 'data/AllSeg.mat'],'AllSeg')

SegLength=sqrt((AllSeg(:,1,1)-AllSeg(:,1,2)).^2 +(AllSeg(:,2,1)-AllSeg(:,2,2)).^2 +(AllSeg(:,3,1)-AllSeg(:,3,2)).^2);
save([TPN 'data/SegLength.mat'],'SegLength')
TotalLength=sum(SegLength);
save([TPN 'data/TotalLength.mat'],'TotalLength')

%% Draw
load([TPN 'D.mat'])
AllSegV=AllSeg/xyum;
AllSegV(:,3,:)=AllSeg(:,3,:)/zum;


Skel=uint8(D)*50;
[ys xs zs] = size(Skel);
Slope = AllSegV(:,:,2)-AllSegV(:,:,1);
Length = sqrt(Slope(:,1).^2 + Slope(:,2).^2 + Slope(:,3).^2);

Divs=2;
Sc=0:Divs;

for l = 1: fix(max(Length))+1
    Get=Length<l;
    Divs=l*10;
    for i = 0: Divs
        Marks=round(AllSegV(Get,:,1)+(Slope(Get,:)/Divs)*i);
        Skel(sub2ind([ys xs zs],Marks(:,1),Marks(:,2),Marks(:,3)))=1000;
    end
end

image(max(Skel,[],3)),pause(.01)
imwriteNp(TPN,Skel,'Skel')

TotalLength

%% NOTES
%{



%}

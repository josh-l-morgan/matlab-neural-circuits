%%identify lengths of dend and count puncta on them

%1 Find Start points
%2 Pick on and start climbing while counting dots. 
%3 When 5(?) um is reached Clip and start over

%4 Insure buffer region exists on either side of segment so that cutting
%segment drop out cant influence density. 



clear all

KPN=GetMyDir
Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;
Buf=1;


%load('cmap.mat')
%set(gcf,'Colormap',cmap)

for k = 1:1%size(Kdir,1)
    clear Branch 
    TPN = [KPN '\' Kdir(k).name '\']
   if exist([TPN 'Use.mat'])
       load([TPN 'Use.mat'])
       load([TPN 'data\AllSegCut.mat'])
       
       NN=Use.NN;
       Mids=Use.Mids;
       Length=Use.Length;
       DotsOnDend=zeros(size(Mids,1),1);
       for d = 1:size(NN,1) %Count Dots on each midpoint
           Dists=dist(Mids,NN(d,:));
           [ID]=find(Dists==0);
           DotsOnDend(ID(1))=DotsOnDend(ID(1))+1;
       end
       
       BranchMap=zeros(size(Mids,1),1); %Map of branch number corresponding to Mids
       BranchNum=1; %Keep track of branch numbers
       
       %% Find Cut length (scaled to mean density)
       load([TPN 'data\Results.mat'])
        meanDD=0;
        for i=1:size(Results.Arbor,2)
            meanDD=meanDD+Results.Arbor(i).DotDend/size(Results.Arbor,2);
        end
        CutBranch=1/meanDD;

            
       Cut=ones(size(Mids,1),1,2)>0;
       
       
%% Clip ends

%% Find all end points
%%Find Connected
EndP=zeros(size(AllSegCut,1),2,'uint8');
for i = 1:size(AllSegCut,1)
    for s=1:2
        Dists=dist(AllSegCut(:,:,:),AllSegCut(i,:,s)); %find distance to all others
        [ID Sing nP] = find3(~Dists.*Cut);
        Con=size(ID,1); %is this Point an end point?
        EndP(ID,nP)=Con==1;
    end
end
[EndID EndnP]=find(EndP);

%% Run Endpoints to clear up to Buf distance


for i = 1:size(EndID,1)     %% Run All Segments

            BranchMapTemp=BranchMap*0;
            %%Initiate Segment
            CheckNext=AllSegCut(EndID(i),:,EndnP(i)); %define segments to begin checking
            BranchLength=0; %initiate Branch length
            DotCount=0; %Start counting Dots
            SumPos=[0 0 0]; %sum positions for final averaging
            MidCount=0; %count mid points for final averaging

            while exist('CheckNext','var') %While there are more segments to check
                CheckNow=CheckNext; %update CheckNow
                clear CheckNext;

                for cn=1:size(CheckNow,1)  %Run all identified segments
                    Dists=dist(AllSegCut(:,:,:),CheckNow(cn,:)); %look for segments with point
                    [ID Sing nP] = find3(~Dists.*Cut);
                    P=~(nP-1)+1; %identify partner index
                    c=0; %Check counter
                    for p=1:size(ID,1) %Add each new point (each branch to list)
                        SumPos=SumPos+Mids(ID(p),:);
                        MidCount=MidCount+1;
                        CheckNext(c+p,:)=AllSegCut(ID(p),:,P(p)); %next segments to check
                        BranchLength=BranchLength+Length(ID(p)); %record length
                        Cut(ID(p),:,:)=0; %Remove segment from future concideration
                        DotCount=DotCount+DotsOnDend(ID(p));
                        BranchMapTemp(ID(p))=BranchNum;
                    end %search all new points
                    c=c+size(ID,1); %extend check list

                end %Run all CheckNows

                if BranchLength>=Buf, break, end %If long enough
            end %Run while new CheckNexts have been found

            %% Record info

end %Run all endpoints

  
       
%% Find Connected
       for rep = 1:3
       for s=1:2
       for i = 1:size(AllSegCut,1)     %% Run All Segments   
           if Cut(i,1,s) %if still present to search
               %Search for Start Point
                Dists=dist(AllSegCut(:,:,:),AllSegCut(i,:,s)); %find distance to all others
                [ID Sing nP] = find3(~Dists.*Cut);
                Con=size(ID,1); %is this Point an end point?
                
            if Con==1 && Cut(i,1,s) %If end point, run segment
                BranchMapTemp=BranchMap*0;               
                %%Initiate Segment
                CheckNext=AllSegCut(i,:,s); %define segments to begin checking
                BranchLength=0; %initiate Branch length
                DotCount=0; %Start counting Dots
                SumPos=[0 0 0]; %sum positions for final averaging
                MidCount=0; %count mid points for final averaging
                
                while exist('CheckNext','var') %While there are more segments to check
                    CheckNow=CheckNext; %update CheckNow 
                    clear CheckNext;
                
                    for cn=1:size(CheckNow,1)  %Run all identified segments
                        Dists=dist(AllSegCut(:,:,:),CheckNow(cn,:)); %look for segments with point
                        [ID Sing nP] = find3(~Dists.*Cut);
                        P=~(nP-1)+1; %identify partner index
                        c=0; %Check counter
                        for p=1:size(ID,1) %Add each new point (each branch to list)
                            SumPos=SumPos+Mids(ID(p),:);
                            MidCount=MidCount+1;
                            CheckNext(c+p,:)=AllSegCut(ID(p),:,P(p)); %next segments to check
                            BranchLength=BranchLength+Length(ID(p)); %record length
                            Cut(ID(p),:,:)=0; %Remove segment from future concideration
                            DotCount=DotCount+DotsOnDend(ID(p));
                            BranchMapTemp(ID(p))=BranchNum;
                        end %search all new points
                        c=c+size(ID,1); %extend check list

                    end %Run all CheckNows
                    
                    if BranchLength>=CutBranch, break, end %If long enough
                end %Run while new CheckNexts have been found
                
                %% Record info
                
                if BranchLength>=CutBranch %If long enough
                    Branch.Pos(BranchNum,:)=SumPos/MidCount;
                    Branch.Length(BranchNum,1)=BranchLength;
                    Branch.DotCount(BranchNum,1)=DotCount;
                    BranchMap=BranchMap+BranchMapTemp;
                    BranchNum=BranchNum+1; %Incriment Branch number
                end
            end %Run point if an end point
            
            
         
            %if mod(i,100)==0, PercentDone=i/(2*size(AllSegCut,1))*100,end
       end %if not cut
       end %Run all segments
       end %Run each side of segment
       end %Repeat to clean up stray bits
          Branch.DD=Branch.DotCount./Branch.Length;
          Branch.meanDD=mean(Branch.DD);  
          Branch.BranchMap=BranchMap;
          
           %% Draw Branchs
           
           cmap =jet(256);
           cmap(1:1,:)=0;
           colormap(cmap)
           
           DrawSeg=fix(Mids)+1;
           DrawBranch=zeros(max(DrawSeg(:,1)),max(DrawSeg(:,2)));
           UseMid=find(BranchMap);
           for ds=1:size(UseMid,1)
                DrawBranch(DrawSeg(UseMid(ds),1),DrawSeg(UseMid(ds),2))=500*Branch.DD(BranchMap(UseMid(ds)))+2;           
           end
           imwrite(DrawBranch,cmap,[TPN 'images\DrawBranchBuf.tif'],'Compression','none')
           image(DrawBranch) ,pause(.3)
           Branch.Info.meanDD=meanDD;
           Branch.Info.CutBranch=CutBranch;
           save([TPN 'BranchSb.mat'],'Branch')
       
       
       
   end %Run cell if Use exists
end %run all cells
%%Find dot Nearest Neighbors

%1 Pick dot on midpoint
%2 Climb to next dot
%3 Record found dot or seg end point
%4 eliminate dots where endpoint is found too early

clear all

KPN=GetMyDir
Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;
CutBranch=20;

load('.\cmap.mat')
set(gcf,'Colormap',cmap)

for k = 1:size(Kdir,1)
    clear NNl NNL
    TPN = [KPN '\' Kdir(k).name '\']; 
   if exist([TPN 'Use.mat'])
       load([TPN 'Use.mat'])
       load([TPN 'data\AllSegCut.mat'])
       
       NN=Use.NN;
       DPos=Use.DPos;
       Mids=Use.Mids;
       Length=Use.Length;
       DotsOnDend=zeros(size(Mids,1),1,2);
       DoDc=cell(size(Mids,1),1,2);
       for d = 1:size(DPos,1) %Count Dots on each midpoint
           Dists=dist(AllSegCut,DPos(d,:));
           [ID h sp]=find3(Dists==min(Dists(:)));
           DID=cell2mat(DoDc(ID(1),h(1),sp(1)));
           DID=[DID d];
           DoDc(ID(1),h(1),sp(1))={DID};
           DotsOnDend(ID(1),h(1),sp(1))=DotsOnDend(ID(1),h(1),sp(1))+1;
       end
         
     
       
       
       %%Find Connected
       Hits=zeros(size(DPos,1),1); HitsN=Hits;
       for i = 1:size(DPos,1)     %% Run All Segments   
                
                %Keep track of segs already traced
                Cut=ones(size(AllSegCut,1),1,2);
           
                %% Find Closest / initiate branching
                Dists=dist(AllSegCut,DPos(i,:));
                [TargS h TargP]=find3(Dists==min(Dists(:)));
                clear CheckNext PathL minDist On segDist PathL
                CheckNext(1,:)=AllSegCut(TargS(1),:,TargP(1));
                On=1; %activate direction
                PathL=0;
                minDist=100;
                segDist=100;
                
                DrawSearch=zeros(200);
                
                while sum(On)
                for nn=1:size(CheckNext,1) %run all directions
                 
                    %%Run Path 
                    if On(nn)   %% If active path, check for dot
                        
                    CheckNow=CheckNext(nn,:); %update CheckNow 
                    Dists=dist(AllSegCut(:,:,:),CheckNow); %look for segments with point
                        [ID Sing nP] = find3(~Dists.*Cut);
                       
                        %% Check for Dot
                        NearDots=[];
                        for nd=1:size(ID,1) 
                            NearDots=[NearDots cell2mat(DoDc(ID(nd), Sing(nd), nP(nd)))];
                        end
                        NearDots=NearDots(NearDots~=i);  %remove current dot
                        if ~isempty(NearDots), HitsN(i)=HitsN(i)+1; end
                        for nd=1:size(NearDots,1)
                            minDist(nn)=min(minDist(nn),dist(DPos(i,:),DPos(NearDots(nd),:)));
                            segDist(nn)=min(segDist(nn),PathL(nn));
                            %On(nn)=0;
                        end
                        
                    end

                    if On(nn)  %%if path still active find next segs
                                           
                        P=~(nP-1)+1; %identify partner index
                        if isempty(ID), On(nn)=0; else  %% Turn off Check
                            for p=1:size(ID,1) %Add each new point (each branch to list)
                                if p==1 %if nothing was found
                                    CheckNext(nn,:)=AllSegCut(ID(p),:,P(p));
                                    PathL(nn)=PathL(nn)+Length(ID(p)); %add to length
                                    Cut(ID(p),:,:)=0; %Remove segment from future concideration
                                    if (PathL(nn)>CutBranch) & (minDist(nn)<100),On(nn)=0; end
                                else
                                    c=size(CheckNext,1)+1;
                                    CheckNext(c,:)=AllSegCut(ID(p),:,P(p)); %next segments to check
                                    On(c)=1;
                                    PathL(c)=PathL(nn)+Length(ID(p));
                                    Cut(ID(p),:,:)=0; %Remove segment from future concideration
                                    minDist(c)=100;
                                    segDist(c)=100;
                                    if (PathL(c)>CutBranch) & (minDist(nn)<100), On(c)=0; end
                                end % if neccessary to grow NNl
                            end %search all new points
                        end %if new point was found

                    end %Run all CheckNows
                    
                  end %Run while new CheckNexts have been found
                end %While there is still an activated path (On)
                %% Record info
                DorL=max(minDist,segDist);
                NNL.DorL(i)=min(DorL);
                NNL.minDist(i)=min(minDist);
                NNL.maxPath(i)=max(PathL);
         
            if mod(i,10)==0, PercentDone=i/(size(DPos,1))*100,end
       end %Run all dots
     
       
       
            scatter(NNL.maxPath,NNL.DorL)
            scatter(NNL.DorL,NNL.minDist)
          
           %% Draw Branchs
           DrawNN=fix(DPos)+1;
           DrawNNL=zeros(max(DrawNN(:,1)),max(DrawNN(:,2)));
           UseNN=find(NNL.maxPath>=5);
           for ds=1:size(UseNN,2)
                DrawNNL(DrawNN(UseNN(ds),1),DrawNN(UseNN(ds),2))=NNL.minDist(UseNN(ds));           
           end
           imwrite(DrawNNL*50,cmap,[TPN 'images\DrawNNL.tif'],'Compression','none')
           image(DrawNNL*10) ,pause(.3)
           save([TPN 'NNL.mat'],'NNL')
       
       
       
   end %Run cell if Use exists
end %run all cells
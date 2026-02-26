%% loading everything

medTis=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\medRes\Analysis\tis.mat');
highTis=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\tis.mat');
medTis=medTis.tis;
highTis=highTis.tis;
medFV='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\medRes\Analysis\fvLibrary\';
highFV='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\fvLibrary\';
mainFV='Y:\karlsRetina\CellNavLibrary_IxQ\Analysis\fvLibrary\';
medObj=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\medRes\Merge\dsObj.mat');
highObj=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Merge\dsObj.mat');
medObj=medObj.dsObj;
highObj=highObj.dsObj;

loadAll=0;
if loadAll
    allSkels=cell(1,length(vgcCidList));
    for i=1:length(vgcCidList)
        curCid=vgcCidList(i);
        skelFileName=['sm_cid' + string(curCid)+'.mat'];
        curSkel=load([skelDir+skelFileName]);
        allSkels{i}=curSkel;
    end

    for i=1:length(allSkels)
        curSkel=allSkels{i};
        curSWC=nep2swc(curSkel.sm.nep);
        allSkels{i}.swc=curSWC;
        
    end
end

%% testing if the skels look okay.

for i=6%1:length(allSkels)
    curSkel=allSkels{i};
    curSM=curSkel.sm;
    curSWC=curSkel.swc;
    curPred=curSWC.pred(curSWC.arbor2swcID)+1;
    curPredInv=zeros(size(curPred));
    curPredInv(curPred>0)=curSWC.swc2arborID(curPred(curPred>0));
    curSM.pred=curPred;
    allEdges=curSM.arbor.edges;
    uniqueNodes=unique(allEdges(:));
    nodeCounts=histcounts(allEdges,uniqueNodes);
    tipIDs=find(nodeCounts==1);
    forkIDs=find(nodeCounts>2);
    distFromRoot=curSM.skel2skel.linDist(tipIDs,1);
    [srtd srtIdx] = sort(distFromRoot,'descend');
    tipIDsSrtd=tipIDs(srtIdx');
    
%     tipTest=[];
%     curNode=1648;
%     its=0;
%     while curNode>0 & its<5000
%         tipTest=[tipTest;curNode];
%         curNode=curPred(curNode);
%         its=its+1;
%     end
    
    usedNodes=[];
    usedEdges=[];
    branchNodes=cell(length(tipIDsSrtd),1);
    branchSize=zeros(length(tipIDsSrtd),1);
    for j=1:length(tipIDsSrtd)
        curTip=tipIDsSrtd(j);
        curBranch=curTip;
        curNode=curTip;
        nextEdge=0;
        stopBool=0;
        counter=0;
        while counter<5000 & curNode>0
            prevEdge=nextEdge(1);
            prevNode=curNode;
            usedNodes=[usedNodes;curNode];
            nextEdge=find(any(allEdges==curNode,2));
            nextEdge=nextEdge(nextEdge~=prevEdge);
            if isempty(nextEdge)
                break
            end
            edgeNodes=allEdges(nextEdge(1),:);
            curNode=edgeNodes(edgeNodes~=prevNode);
            curBranch=[curBranch;curNode];
            %if ismember(curNode,usedNodes)
            %    break
            %end
            counter=counter+1;
        end
        branchNodes{j}=curBranch;
        branchSize(j)=length(curBranch);
    end
        
    if 1
        figures{i}=figure();
        figStruct=compareMorph(figures{i},curSM.cid,highFV);
        title(curSM.cid)
        hold on
        for k=1:length(branchNodes)
            curBranch=branchNodes{k};
            branchPos=curSM.nep.pos(curBranch,:);
            figStruct.skel(k)=plot3(branchPos(:,3),branchPos(:,1),branchPos(:,2),'c-','LineWidth',1);
            %credge=allEdges(k,:);
            %plot3(curSM.nep.pos(credge,1),curSM.nep.pos(credge,2),curSM.nep.pos(credge,3));
            %drawnow();
        end
        figStruct.patches(1).FaceAlpha=0.1;
        figStruct.patches(1).FaceColor=[.7 0 .7];
        %figStruct.skel.LineWidth=3;
    end
end


%% testing the bwskel function
curVox=getCidVox(14,1,highObj,highTis);
curVox=curVox{1};
curVox=curVox(curVox(:,1)>1,:);
bbox=[min(curVox(:,1)),max(curVox(:,1));min(curVox(:,2)),max(curVox(:,2));min(curVox(:,3)),max(curVox(:,3))];
boxDims=bbox(:,2)-bbox(:,1)+1;
miniVox=curVox-(bbox(:,1)'-1);
volImg=zeros(boxDims(1),boxDims(2),boxDims(3),'logical');
indList=sub2ind(size(volImg),miniVox(:,1),miniVox(:,2),miniVox(:,3));
volImg(indList)=1;
hsize=9;sigma=.5;
h = fspecial3('gaussian',hsize,sigma);
%testImg=double(volImg);
volImgBlurg=convn(volImg,h,'same');
volImgThrsh=volImgBlurg>0.1;


if 0
testFig=figure();
hold on
for j=size(volImg,3):-3:1
    image=zeros(size(volImg,1),size(volImg,2),3);
    image(:,:,[1 3])=image(:,:,[1 3])+volImg(:,:,j).*0.5;
    image(:,:,[2 3])=image(:,:,[2 3])+volImgBlurg(:,:,j).*0.5;
    image(:,:,[1 2])=image(:,:,[1 2])+volImgThrsh(:,:,j).*0.5;
    imshow(image);
    drawnow();
    pause(0.01);
end
testFig2=figure();
hold on
end

%volImgThrsh=volImgBlurg>0.3;
%[fac,vert]=isosurface(volImg);
bwsk=bwskel(volImg);
bwskThrsh=bwskel(volImgThrsh);

skelThrshSurf=isosurface(bwskThrsh);
[x,y,z]=ind2sub(size(bwskThrsh),find(bwskThrsh>0));
bwskThrshCoords=horzcat(x,y,z);
figure(); scatter3(x,y,z);

figure(); p=patch(skelThrshSurf); p.FaceAlpha=0.25; p.FaceColor=[1 0 1];

branchPts=bwmorph3(bwskThrsh,'branchpoints');
endPts=bwmorph3(bwskThrsh,'endpoints');

for j=size(bwsk,3):-1:1
    imshow(bwsk(:,:,j));
    drawnow();
    pause(0.01);
end
skelfv=struct();
skelfv.Vertices=bvert;
skelfv.Faces=bfac;
skelpat=patch(skelfv);




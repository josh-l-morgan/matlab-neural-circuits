%fancifig
rgcList=type2cid({'rgc'},{'4ow'},curTis);
alphList=rgcList{1};

alphFig=figure();
hold on
alphPat=cell(length(alphList),1);
for i=1:length(alphList)
    curAlph=alphList(i);
    curFV=load([fvDir num2str(curAlph) '.mat']);
    curPat=patch(curFV.fv);
    curPat.FaceColor='flat';
    curPat.EdgeColor='none';
    alphPat{i}=curPat;
end
title('patch each cell')

%% Plot as one thing.
% Karl, finish this later. You need to add the size(1) of the object to the 
% face IDs or it will not be able to draw the right faces. 
alphFig2=figure();
hold on
superFV.fv=struct();
superFV.fv.vertices=[0 0 0];
superFV.fv.faces=[0 0 0];
facePad=0;
for i=1:length(alphList)
    curAlph=alphList(i);
    curFV=load([fvDir num2str(curAlph) '.mat']);
    superFV.fv.vertices=vertcat(superFV.fv.vertices,curFV.fv.vertices);
    facePad=facePad+size(curFV.fv.vertices,1);
    superFV.fv.faces=vertcat(superFV.fv.faces,(curFV.fv.faces+facePad));
end
superFV.fv.vertices(1,:)=[];
superFV.fv.faces(1,:)=[];
superPat=patch(superFV.fv);
superPat.EdgeAlpha=0;
superPat.EdgeColor=[0 0 0];
title('single patch')
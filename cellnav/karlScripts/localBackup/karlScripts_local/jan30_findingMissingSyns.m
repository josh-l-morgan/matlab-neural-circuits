oldTis=load('C:\work\liteLib\tis.mat');
highTis=load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\tis.mat');
oldTis=oldTis.tis;
highTis=highTis.tis;

synCheck=repmat(9,length(oldTis.syn.pos(:,1)),11);
for synIt=1:length(oldTis.syn.pos(:,1))
    curPos=oldTis.syn.pos(synIt,:);
    curEdges=oldTis.syn.edges(synIt,:);
    mat=find(pdist2(curPos,highTis.syn.pos)<.2);
    if mat
    matEdges=highTis.syn.edges(mat,:);
    matPos=highTis.syn.pos(mat,:);
    if sum(curEdges([1 2])==matEdges(:,[1 2]),'all')>0
        goodRow=find(sum(curEdges([1 2])==matEdges(:,[1 2]),2)==max(sum(curEdges([1 2])==matEdges(:,[1 2]),2)));
        synCheck(synIt,:)=[curPos matPos(goodRow(1),:) curEdges([1 2]) matEdges(goodRow(1),[1 2]) 1];
    else
        synCheck(synIt,:)=[curPos matPos curEdges([1 2]) matEdges([1 2]) 2];
    end
    else
        synCheck(synIt,:)=[curPos 0 0 0 curEdges([1 2]) 0 0 3];
    end
end


length(find(oldTis.syn.edges(:,1)==2002))
length(find(highTis.syn.edges(:,1)==2002))

tab1=tabulate(oldTis.syn.edges(oldTis.syn.edges(:,1)==2002,2));
tab2=tabulate(highTis.syn.edges(highTis.syn.edges(:,1)==2002,2));

[~,ord1]=sort(tab1(:,2),'descend');
[~,ord2]=sort(tab2(:,2),'descend');

searchCid=5;
locA=oldTis.syn.pos(oldTis.syn.edges(:,2)==searchCid&oldTis.syn.edges(:,1)==2002,:);
locB=uint16(locA(:,[2 1 3]).*[250 250 25]);

tab1srtd=tab1(ord1,:);
tab2srtd=tab2(ord2,:);



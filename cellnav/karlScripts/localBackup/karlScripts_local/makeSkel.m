function sklout=makeSkel(cid,dsObj,curTis)

curVox=getCidVox(cid,1,dsObj,curTis);
curVox=curVox{1};
curVox=curVox(curVox(:,1)>1,:);
bbox=[min(curVox(:,1)),max(curVox(:,1));min(curVox(:,2)),max(curVox(:,2));min(curVox(:,3)),max(curVox(:,3))];
boxDims=bbox(:,2)-bbox(:,1)+1;
miniVox=curVox-(bbox(:,1)'-1);
volImg=zeros(boxDims(1),boxDims(2),boxDims(3),'logical');
indList=sub2ind(size(volImg),miniVox(:,1),miniVox(:,2),miniVox(:,3));
volImg(indList)=1;

volImgDil=volImg;
dilator=5;
for q=ceil(dilator/2):size(volImg,3)-floor(dilator/2)
    wind=volImg(:,:,[q-floor(dilator/2):q+floor(dilator/2)]);
    %windMid=volImg(:,:,q);
    %windMode=mode(wind,3);
    %windMid(windMode)=1;
    volImgDil(:,:,q)=max(wind,[],3);
end

testVols={volImg,volImgDil};

if 0
testFig=figure();
hold on
for j=size(volImg,3):-3:1
    image=zeros(size(volImg,1),size(volImg,2),3);
    image(:,:,[1 3])=image(:,:,[1 3])+volImg(:,:,j).*0.5;
    image(:,:,[2 3])=image(:,:,[2 3])+volImgDil(:,:,j).*0.5;
    %image(:,:,[1 2])=image(:,:,[1 2])+volImgThrsh(:,:,j).*0.5;
    imshow(image);
    drawnow();
    pause(0.01);
end
testFig2=figure();
hold on
end

testSkels={};
testBPs={};
testTips={};
for u=1:length(testVols)
    curVol=testVols{u};
    curSk=bwskel(curVol);
    testSkels{u}=curSk;
    curSkBPs=bwmorph3(curSk,'branchpoints');
    curSkTPs=bwmorph3(curSk,'endpoints');
    testBPs{u}=curSkBPs;
    testTips{u}=curSkTPs;
    [x,y,z]=ind2sub(size(curSk),find(curSk==1));
    testSkelVox{u}=horzcat(x,y,z);
end

%% plot
figure();
scatter3(testSkelVox{1}(:,1),testSkelVox{1}(:,2),testSkelVox{1}(:,3),5,'m.');
testPt=[ 197   883   401 ];


end


%% other 
%figure();
%hold on
skelPos={};
bpPos={};
tipPos={};
for o=1:length(testSkels)
    curSk=testSkels{o};
    curBP=testBPs{o};
    curTip=testTips{o};
    [x,y,z]=ind2sub(size(curSk),find(curSk>0));
    [xBP,yBP,zBP]=ind2sub(size(curBP),find(curBP>0));
    [xT,yT,zT]=ind2sub(size(curTip),find(curTip>0));
    skelPos{o}=horzcat(x,y,z);
    bpPos{o}=horzcat(xBP,yBP,zBP);
    tipPos{o}=horzcat(xT,yT,zT);
    %scatter3(x,y,z,2,'.');
end


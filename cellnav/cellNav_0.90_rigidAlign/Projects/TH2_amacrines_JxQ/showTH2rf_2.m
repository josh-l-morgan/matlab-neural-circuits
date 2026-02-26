

if 0
    global glob tis
    SPN =  [glob.datDir 'Analysis\Data\'];

    volName = glob.vol.activeName;
    smDir = [glob.datDir 'Volumes\' volName '\Analysis\SMs\']


    showCell = [10];
    cid = showCell(1);
    smName = sprintf('%ssm_cid%d.mat',smDir,cid)
    load(smName)
end


lc = 100;
bipRFsigma = 50;
bipGau = fspecial('gaussian',bipRFsigma*8,bipRFsigma);
bipGau = bipGau * 1/max(bipGau(:));

yOffset = 25;
xOffset = 125;


nep = sm.nep;
pos = nep.pos;

minPos = min(pos,[],1);
maxPos = max(pos,[],1);
fSize = round(maxPos + 50);
field = zeros(fSize(1),fSize(2));

d = sm.syn2Skel.syn2SkelDist;

W = exp(-d/lc); % Apply length constant


%% Assign RF to syn
syn = sm.syn;
isRib = sm.syn.synType == 2;
pre = sm.syn.pre;

usePre = find(syn.synType == 2);
fStack = zeros(fSize(1),fSize(2),length(usePre));

for s = 1:length(usePre)
    synID = usePre(s);
    
    preCid = pre(usePre(s));
    
    %targ = find(
    


    sField = field;
    sField(round(syn.pos(synID,1)),round(syn.pos(synID,2))) = 1;
    fField = imfilter(sField,bipGau);
    fStack(:,:,s) = fField;
    if 1
    clf
    image((fField+sField)*200)
    hold on
    scatter(pos(:,2),pos(:,1),'.','k')
    scatter(syn.pos(synID,2),syn.pos(synID,1),'r','filled')
    hold off
    drawnow
    end
end

rfSum = sum(fStack,3);
image(rfSum*4)

%% run
useN = 1:10:size(d,2);
yMax = zeros(length(useN),1);
xMax = yMax;
for uN = 1:length(useN)
    n = useN(uN);
    disp(sprintf('%d of %d',n,size(d,2)))
    clear sW
    nW(1,1,:) = W(usePre,n);
    nStack = fStack .* repmat(nW,[fSize(1) fSize(2) 1]);
    nSum = sum(nStack,3);

    [y x] = find(nSum == max(nSum(:)),1);
    yMax(uN) = y;
    xMax(uN) = x;

    if ~mod(uN-1,round(length(useN)/10))
        clf
        image(nSum * 50/mean(nSum(:)))
        hold on
        scatter(pos(:,2),pos(:,1),'.','k')
        scatter(pos(n,2),pos(n,1),100,'g','filled')
        scatter(x,y,100,'r','filled')
        hold off
        drawnow
    end
    %         yMax(uN) = pos(n,1);
    %         xMax(uN) = pos(n,2);
end


col = jet(100);

yInd = round((yMax-yOffset)/200 * 100);
yInd(yInd<1) = 1;
yInd(yInd>100) = 100;

xInd = round((xMax-xOffset)/200 * 100);
xInd(xInd<1) = 1;
xInd(xInd>100) = 100;

%% Draw color bars

clf


subplot(2,1,1)

scatA = scatter3(pos(:,1),pos(:,2),pos(:,3),3,'k','o','filled',...
    'markerfacealpha',.2,'markerfacecolor','flat');
hold on
scat1 = scatter3(pos(useN,1),pos(useN,2),pos(useN,3),'r','o','filled');
view(0,90)
scat1.CData = col(yInd,:);
hold on
yBar = [1:200]+yOffset;
scatY = scatter3(yBar,yBar*0+xOffset,yBar*0+maxPos(3),'o','filled');
scatY.CData = col(round((yBar-yOffset)/2),:);
hold on

subplot(2,1,2)

scatB = scatter3(pos(:,1),pos(:,2),pos(:,3),3,'k','o','filled',...
    'markerfacealpha',.2,'markerfacecolor','flat');
hold on
scat2 = scatter3(pos(useN,1),pos(useN,2),pos(useN,3),'r','o','filled');
view(0,90)
scat2.CData = col(xInd,:);
xBar = [1:200]+xOffset;
scatX = scatter3(xBar*0+yOffset,xBar,xBar*0+maxPos(3),'o','filled');
scatX.CData = col(round((xBar-xOffset)/2),:);
hold off


%% Plot shift vectors
clf
scatA = scatter3(pos(:,1),pos(:,2),pos(:,3),3,'k','o','filled',...
    'markerfacealpha',.2,'markerfacecolor','flat');
hold on
axis 'equal'
set(gca,'clipping', 'off')
view(0,90)

for i = 1:length(useN)
    
    startP = pos(useN(i),:);
    stopP = [yMax(i) xMax(i) pos(useN(i),3)];
    distP = sqrt((startP(1)-stopP(1)).^2 + (startP(2)-stopP(2)).^2 + ...
        (startP(3)-stopP(3)).^2);
    distInd = round(distP);
    distInd(distInd<1) = 1;
    distInd(distInd>100) = 100;
    
    plot3([startP(1) stopP(1)],[startP(2) stopP(2)],[startP(3) stopP(3)],...
        'linewidth',3,'color',col(distInd,:));
    scatter3(startP(1),startP(2),startP(3),'k','o','filled')
    %scatter3(stopP(1),stopP(2),stopP(3),'r','o','filled')



end
hold off



%% Save eps
if 0

    f = gcf;

    pause(.01)
    saveDir = uigetdir;
    imageName = [saveDir '\y100umLCb.eps'];
    apos= get(f,'Position');
    set(f,'PaperUnits','points','PaperPosition', apos);
    set(f, 'InvertHardCopy', 'off');
    print(gcf,[imageName],'-depsc','-r256','-opengl');
    print(gcf,[imageName],'-depsc','-r256','-opengl','-noui', '-vector');
    pause(.03)




end




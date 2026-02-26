
%% colect data on cells

cellSynCount = zeros(length(s),length(useVgc));
clear synPos cellSynPos
countG = zeros(length(s),1);
for i = 1:length(s)
    synPos{i} = [];
end
allPos = [];
allVoxPos = [];
allLengths = [];
vArborLengths = [];
allDepths = [];
% T = table;
% T.Properties.VariableNames = {'cid' 'Area' 'MajorAxisLength' ...
%     'MinorAxisLength' 'Orientation'};
clear Area MajorAxisLength MinorAxisLength Orientation
for v = 1:length(useVgc)

    cid = useVgc(v);
    useSM(v) = 1;
    disp(sprintf('loading sm %d, %d of %d',cid,v,length(useVgc)))
    sm = sms(v).sm;%load([smDir fileName],'syn','syn2Skel','nep');

    syn = sm.syn; % get syn information
    synD = sm.syn2Skel.syn2SynDist;
    skelD = sm.syn2Skel.syn2SkelDist;
    nodeLengths = sm.nep.props.nodeLength;  


    allVoxPos = cat(1,allVoxPos,sm.arbor.subs);
    allPos = cat(1,allPos,sm.nep.pos);
    allLengths = cat(1,allLengths,nodeLengths);

      

    %%Lengths

    cPos = sm.nep.pos;
    [zGCL,zINL,nDepth1] = getIPLdepth(cPos(:,3),cPos(:,1),cPos(:,2),[],[]);
    notCB = nDepth1 > 0.1;
    vArborLengths(v) = sum(nodeLengths(notCB));
    allDepths = cat(1,allDepths,nDepth1);

    
    ys = ceil(max(cPos(:,1)))+100;
    xs = ceil(max(cPos(:,2)))+100;
    cI = zeros(ys,xs);
    cI(sub2ind([ys xs],round(cPos(:,1)),round(cPos(:,2)))) = 1;
    image(cI*1000)
    
    cProps = regionprops(cI,'Area','Circularity','ConvexArea',...
        'Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength',...
        'Orientation');
    
    Area(v,1) = cProps.ConvexArea;
    MajorAxisLength(v,1) = cProps.MajorAxisLength;
    MinorAxisLength(v,1) = cProps.MinorAxisLength;
    Orientation(v,1) = cProps.Orientation;  

end

Area = sprintf('%0.0f\n',Area)
MajorAxisLength = sprintf('%0.0f\n',MajorAxisLength)
MinorAxisLength = sprintf('%0.0f\n',MinorAxisLength)
Orientation = sprintf('%0.0f\n',Orientation)


%% Get depth data

depthBins = 0:.001:1;
binW = .025;
clear lengthAtDepth  atDepth
for i = 1:length(depthBins)-1
    b1 = depthBins(i)-binW;
    b2 = depthBins(i+1)+binW;
    isDeep = (allDepths>b1) & (allDepths<=b2);
    lengthAtDepth(i) = sum(allLengths(isDeep));
    atDepth(i) = mean([b1 b2]);
end

plot(atDepth,lengthAtDepth)




%%Find upper and lower bounds
bound1 = 0.025;
bound2 = 0.975;
lowBound = [];
highBound = [];
totLength = sum(lengthAtDepth);
trackLength = 0;
for i = 1:length(lengthAtDepth)
    trackLength = trackLength + lengthAtDepth(i);
    if isempty(lowBound)
        if trackLength/totLength > bound1
            lowBound = depthBins(i);
        end
    end
    if isempty(highBound)
        if trackLength/totLength > bound2
            highBound = depthBins(i);
        end
    end

end

%T = table(useVgc,Area(1:end-1),MajorAxisLength(1:end-1),MinorAxisLength(1:end-1),Orientation(1:end-1))







%Generate violyn scatter plots of nearby synapses. Grouped and colored by type.

global glob tis

SPN = [glob.datDir 'Analysis\Data\preproc\'];
load([SPN 'COI.mat']);
load([SPN 'NOI.mat']);

smDir = [glob.dir.Volumes  glob.vol.activeName '\Analysis\SMs\'];

useVgc = COI.vgcCids;%[2 3 4 14];

distBin = 1;
synGroups = [1 2 3 4 5 0];
synTypes = [1 2 3 4 5 0];

bRange = [0 : distBin : 50];
g2g = zeros(length(synGroups),length(synTypes),length(bRange)-1);
nearLength = g2g;
countG = zeros(length(synGroups),1);
allSynTypes=[];
for v = 1:length(useVgc)
    
    cid = useVgc(v);
    
    fileName = sprintf('smx_cid%d.mat',cid);
    useSM(v) = 1;
    disp(sprintf('loading sm %d, %d of %d',cid,v,length(useVgc)))
    sm = load([smDir fileName],'syn','syn2Skel','nep');
    
    syn = sm.syn; % get syn information
    synD = sm.syn2Skel.syn2SynDist;
    skelD = sm.syn2Skel.syn2SkelDist;
    synD(sub2ind(size(synD),1:size(synD),1:size(synD))) = inf;
    nodeLengths = sm.nep.props.nodeLength;
    
    %synformation parsing for subtype
    preCheck=cid2type(syn.pre,curTis);
    syn.preSub=preCheck{3}';
    postCheck=cid2type(syn.post,curTis);
    syn.postSub=postCheck{3}';
    
    %% sort types
    synType = zeros(length(syn.preClass),1);
    %synType( (syn.preClass == 7) | (syn.synType == 2) ) = 1; % bip
    %synType( (syn.synType == 1) & (syn.post == cid) ) = 1; % amc in
    %synType( (syn.preClass == 0) & (syn.post == cid) ) = 1; % amc in
    synType( (syn.post == cid) & (syn.preClass == 7) & ismember(syn.preSub,bpcOFFsubs)) = 1;
    synType( (syn.post == cid) & (syn.preClass == 7) & ismember(syn.preSub,bpcONsubs)) = 2;
    synType( (syn.post == cid) & (syn.preClass == 7) & ismember(syn.preSub,bpcOFFsubs)) = 1;
    synType( (syn.post == cid) & (syn.preClass == 7) & ismember(syn.preSub,bpcONsubs)) = 2;
    synType( (syn.post == cid) & (syn.preClass == 7) & ismember(syn.preSub,bpcOFFsubs)) = 1;
    synType( (syn.post == cid) & (syn.preClass == 7) & ismember(syn.preSub,bpcONsubs)) = 2;
    %synType( (syn.pre == cid) & (syn.postClass == 1) ) = 3; % rgc
    %synType( (syn.pre == cid) & (syn.postClass == 8) ) = 4; % amc out
    synType( (syn.pre == cid) & (syn.postClass == 1) & (syn.postSub == 21)) = 3; %51
    synType( (syn.pre == cid) & (syn.postClass == 1) & (syn.postSub == 23)) = 4; %4i
    synType( (syn.pre == cid) & (syn.postClass == 1) & (syn.postSub == 24)) = 5; %4ow
    
    allSynTypes=vertcat(allSynTypes,synType);
    %synType( (syn.pre == cid) & (syn.postClass == 1) & (syn.postSub == 21)) = 5; %
    
    %%
    
    for g = 1:length(synGroups);
        disp(sprintf('checking group %d, %d of %d',synGroups(g),g,length(synGroups)))
        sG = find(synType == synGroups(g));
        countG(g) = countG(g) + length(sG);
        synDG = synD(sG,:);
        skelDG = skelD(sG,:);
        for b = 1:(length(bRange)-1)
            
            [ y x ] = find(((synDG>=bRange(b)) & (synDG<bRange(b+1))));
            if ~isempty(x)
                foundTypes = synType(x);
                for hT = 1:length(synTypes);
                    g2g(g,hT,b) =  g2g(g,hT,b) + sum(foundTypes==synTypes(hT));
                end
            end
            
            [ y x ] = find(((skelDG>=bRange(b)) & (skelDG<bRange(b+1))));
            nearLength(g,:,b) = nearLength(g,:,b) + sum(nodeLengths(x));
            
            
        end
    end
end

%% Scale my number of gs found
g2g = g2g ./ repmat(countG,[1 size(g2g,2) size(g2g,3)]);
nearLength = nearLength ./ repmat(countG,[1 size(g2g,2)]);

g2gN = g2g ./ nearLength;


%% Display group to group

clf
hold on
tCol = [1 0 0; 0 1 0 ; 0 0 1; 1 1 0; 1 0 1; 0 1 1]*.7;
sumT = squeeze(sum(g2g,2));
maxT = max(sumT(:));

titleMat={'OFFbpc','ONbpc','4i','4ow','W3','other'};

for g = 1:3 %1:size(g2g,1)
    for t = 1:2 %1:size(g2g,2)
        %subplot(size(g2gN,1),size(g2gN,2),(g-1)*size(g2gN,1)+t)
        subplot(3,2,(g-1)*2+t)
        plot(squeeze(g2gN(g+2,t,:)))
        ylim([0 .1])
        xlim([0 10])
        title([titleMat{t} ' x ' titleMat{g+2}]);
    end
end

pause(3)
%% Scatter syn
clf
hold on

for t = 1:size(g2gN,2)
    scatStore(t).X = [];
    scatStore(t).Y = [];
end

g2gN2 = g2gN * 100;

spacer = .6/max(g2gN2(:));
for g =  1:size(g2gN2,1)
    
    for b = 1:min(10,size(g2gN2,3))
        tVec = g2gN2(g,:,b);
        sumVec = sum(tVec);
        for t = 1:length(tVec)
            if tVec(t)
                n = tVec(t);
                startX = (n*spacer)/2 * -1;
                shiftX = [startX:spacer:startX+(n * spacer)] + g;
                shiftY = (shiftX*0+bRange(b)) + ((t/length(tVec))*distBin*.6);
                %scatter(shiftX,bRange(b),10,tCol(t,:))
                scatStore(t).X = cat(2,scatStore(t).X,shiftX);
                scatStore(t).Y = cat(2,scatStore(t).Y,shiftY);
                
            end
        end
        
    end
end
for t = 1:size(g2gN2,2)
    
    scatter(scatStore(t).X,scatStore(t).Y,30,tCol(t,:),'filled')
    hold on
    
end


%% RGC to BPC
clf
plot(squeeze(g2gN(3,1,:)),'g')
hold on
plot(squeeze(g2gN(3,3,:)),'b')
hold off
ylim([0 .15])

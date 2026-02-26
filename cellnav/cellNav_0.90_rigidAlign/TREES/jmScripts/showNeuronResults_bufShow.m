
fig = figure;
runCid = 4;
recordAllNodes = 1;
storeAllOut = 0;
w = warning ('on','all')

%% Set up directories
expName = 'standardVal01';
nnDir = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\nn\';
expDir = [nnDir expName '\'];
cellDir = sprintf('%scid%d\\',expDir,runCid);
expInfoName = sprintf('%sexpInfo.mat',cellDir);


if ~exist(nnDir,'dir'),mkdir(nnDir); end
if ~exist(expDir,'dir'),mkdir(expDir); end
if ~exist(cellDir,'dir'),mkdir(cellDir); end


%% Get skeleton
swcDir = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\swc\';

swcFile = sprintf('%scid%d.swc',swcDir,runCid);
mtrFile = sprintf('%scid%d.mtr',swcDir,runCid);


%% Get nep
nepFile = sprintf('%ssm2nrn_cid%d.mat',swcDir,runCid);
load(nepFile);
%load('Z:\Active\morganLab\MATLAB\TREES\morphos\sm2nrn_cid20.mat')


%% Read in the result
load(expInfoName);
expInfo
numExp = expInfo.numExp;

c = 0;

dCell = dir([cellDir 'res*.mat']);
rNames = {dCell.name};
dLength = length(rNames);
load([cellDir rNames{1}])
synIDs = nnRes.synIDs;
numSyn = length(synIDs);
numTime = size(nnRes.synV,2);
numNode = length(nnRes.maxV);

synVs = zeros(dLength,numSyn,numTime);
maxVs = zeros(dLength,numNode);

for e = 1:numExp
    expFileName = sprintf('%sres%d.mat',cellDir,e);
    if exist(expFileName,'file')
        c = c+1;
        sprintf('found exp %d',e)
        load(expFileName)
        synVs(c,:,:) =  nnRes.synV;
        maxVs(e,:) = nnRes.maxV;
    end
end



%% plot the result (Vmem at soma and synapse and synaptic current)
clf
subplot(2,1,1)
hold all
for s = 1:length(synIDs)
    plot(out.t,out.record{1}.cell.v{synIDs(s)})
end
legend('Synapse')
ylabel('Membrane potential [mV]')
ylim([-85,0])
xlim([0,50])
subplot(2,1,2)
plot(out.t,out.record{1}.AlphaSynapse.i{synIDs(end)})    % plot time vs synaptic current
xlim([0,50])
ylabel('Synaptic current [nA]')
xlabel('Time [ms]')


%%
synVs(e,:,:) = synV(uSynIdx,:);
synMax = max(synVs,[],3);
synDev = synMax;
synDev = synDev - median(synDev(:));
synDev(synDev<0) = 0;
sumDev = sum(synDev,3);
allProp = sumDev;

for i = 1:numSyn

    prop = allProp(i,:);
    prop = prop-min(prop);
    clf
    scatter3(tr.X,tr.Y,tr.Z,'.','k')
    axis 'equal'
    hold on
    synScat = scatter3(synPos(:,1),synPos(:,2),synPos(:,3));
    set(synScat,'CData',prop)
    scatter3(tr.X(nearest),tr.Y(nearest),tr.Z(nearest),'p','g')


end

%%
if recordAllNodes


    props = maxVs;
    props = props-min(props(:));
    props = ceil(props * 255/(max(props(:))));
    cmap = jet(  );
    colormap(cmap)

    for i = 1:numSyn

        clf
        plot_tree(tree{1},props(i,:)');
        hold on
        scatter3(synPos(i,1),synPos(i,2),synPos(i,3),250,'k')

        colorbar   % plot tree with region code
        set(gca,'clipping','off')
        pause


    end
end




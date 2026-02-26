

runCids = [ 2 3 4 5 13 14 20]

for v = 1:length(runCids)

%% Define variables
runCid = runCids(v);
recordAllNodes = 1;
storeAllOut = 0;
w = warning ('on','all')
expName = 'testONOFF_2';

%% Set up directories

nnDir = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\nn\';
expDir = [nnDir expName '\'];
cellDir = sprintf('%scid%d\\',expDir,runCid);
expInfoName = sprintf('%sexpInfo.mat',cellDir);

if ~exist(nnDir,'dir'),mkdir(nnDir); end
if ~exist(expDir,'dir'),mkdir(expDir); end
if ~exist(cellDir,'dir'),mkdir(cellDir); end

%% Make figure
set(0,'DefaultFigureWindowStyle','docked')
fig = figure;


%% Make default neuron

neuron = [];                                                % clear neuron structure
neuron.params.v_init = -60;                                 % starting membrane potential [mV] of all cells
neuron.params.dt = 0.025;                                   % integration time step [ms]
neuron.params.tstop = 50;                                  % stop simulation after this (simulation) time [ms]
neuron.params.prerun = -200; %400                           % add a pre runtime [ms] to let system settle
neuron.params.celsius = 35;                                 % temperature [celsius]
neuron.params.nseg = 'dlambda';                             % the dlambda rule is used to set the amount of segments per section. Alternatively, a fixed number can be entered or a string 'EachX' with X being the interval between segments [micron]
neuron.params.accuracy = 0;                                 % optional argument if number of segments should be increased in the soma and axon

%% Set neuron model mechanisms

neuron.experiment = 'test';                                                 % give your simulation run a name (necessary when using advanced t2n protocols)
for t = 1                                                                   % preconfiguration to loop through several morphologies, now only one morphology is considered
    neuron.mech{t}.all.pas = struct('g',0.0003,'Ra',100,'e',-60,'cm',1);    % add passive channel to all regions and define membrane capacity [µF/cm²], cytoplasmic resistivity [Ohm cm] and e_leak [mV]
    neuron.mech{t}.soma.pas = struct('g',0.0004,'Ra',100,'e',-60,'cm',1);    % add passive channel to somatic regions (will overwrite the "all" definition) and define membrane capacity [µF/cm²], cytoplasmic resistivity [Ohm cm] and e_leak [mV]
    neuron.mech{t}.all.k_ion.ek = -90;
    neuron.mech{t}.all.na_ion.ena = 50;
    %neuron.mech{t}.soma.hh = struct('gnabar',0.25,'gkbar',0.036,'gl',0);     % add Hodgkin-Huxley Sodium channel only to soma
    neuron.record{t}.cell = struct('node',1,'record','v');                   % record voltage "v" from first node (i.e. the soma)
end

%% Get skeleton

swcDir = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\swc\';
swcFile = sprintf('%scid%d.swc',swcDir,runCid);
mtrFile = sprintf('%scid%d.mtr',swcDir,runCid);

if exist(mtrFile,'file')
    tree = load_tree(mtrFile);
else
    tree{1} = load_tree(swcFile);
    tree{1}.rnames{1} = 'dendrites';
    tree{1}.rnames{2} = 'subtree';
    tree{1}.Ri = 100;
    tree{1}.Gm = 100 * 10^(-6);
    tree{1}.Cm = 1;
    % tree{1}.R(1:2) = 3;     % make the first two nodes a new region
    % tree{1}.rnames{3} = 'soma';   % name the region "soma"
    % tree{1}.D(1:2) = 10;        % increase diameter of the somatic nodes
    tree = t2n_writeTrees(tree,[],mtrFile);                 % transform tree to NEURON morphology (.hoc file). This only has to be done once for each morphology
end
tree{1}.rnames{1} = 'dendrites';
tree{1}.rnames{2} = 'subtree';
tree{1}.Ri = 100;
tree{1}.Gm = 100 * 10^(-6);
tree{1}.Cm = 1;

%% Fix diameters
minDiam = 0.1;
tree{1}.D(tree{1}.D<minDiam) = minDiam;
if 0
    tree{1}.D = tree{1}.D * 0 + median(tree{1}.D);
end


%% Show cell
clf
plot_tree(tree{1},tree{1}.R);colorbar   % plot tree with region code
axis off

%% Get nep

nepFile = sprintf('%ssm2nrn_cid%d.mat',swcDir,runCid);
load(nepFile);

syn = sm2nrn.syn;
numNodes = length(tree{1}.X);
synPos = syn.pos;
numSyn = size(synPos,1);
tr = tree{1};
minDists = zeros(numSyn,1);
nearest = zeros(numSyn,1);
for s = 1:numSyn
    dists = sqrt((tr.X-synPos(s,1)).^2 + (tr.Y-synPos(s,2)).^2 + (tr.Z-synPos(s,3)).^2);
    minDists(s) = min(dists);
    nearest(s) = find(dists == minDists(s),1);
end

scatter3(tr.X,tr.Y,tr.Z,'.','k')
axis 'equal'
hold on
scatter3(synPos(:,1),synPos(:,2),synPos(:,3),'r')
scatter3(tr.X(nearest),tr.Y(nearest),tr.Z(nearest),'p','g')

hold off

%% Get synapse information

SPN = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Analysis\Data\preproc\'
load('Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\tis.mat')
load([SPN 'COI.mat']);

useOnGroup = [6 7 8 9 10 11 12 14 17 20 21 23];
useOffGroup = [1 2 3 4 5 15 18 19];
onBips = [];
offBips = [];
for i = 1:length(tis.cells.type.typeID)
    if tis.cells.type.typeID(i) == 7
        if intersect(useOnGroup,tis.cells.type.subTypeID(i))
            onBips = [onBips tis.cells.cids(i)];
        elseif intersect(useOffGroup,tis.cells.type.subTypeID(i))
            offBips = [offBips tis.cells.cids(i)];
        end
    end
end

ONsyn = [];
OFFsyn = [];
for i = 1:length(syn.pre);
    if intersect(onBips,syn.pre(i))
        ONsyn = [ONsyn i];
    elseif intersect(offBips,syn.pre(i))
        OFFsyn = [OFFsyn i];
    end
end

synIDs = nearest;
activate{1} = synIDs(OFFsyn);
activate{2} = synIDs(ONsyn);

clf
hold on
view(0,0)
s1 = scatter3(tr.X,tr.Y,tr.Z,1,'.','markeredgecolor',[.6 .6 .6]);
s2 = scatter3(tr.X(activate{1}),tr.Y(activate{1}),tr.Z(activate{1}),50,'o','markeredgecolor',[.8 0 0]);
s3 = scatter3(tr.X(activate{2}),tr.Y(activate{2}),tr.Z(activate{2}),50,'o','markeredgecolor',[0 .8 0]);



%% Make synapse experiment

nneuron = neuron;                                                         % copy standard neuron structure
plen = Pvec_tree(tree{1});                                                % get path length to soma at each node of morphology
%[~,synIDs] = max(plen);                                                   % search for the most far away point in the morpholgy

%%Force unique
[uSynIDs a uSynIdx] = unique(synIDs);

%%Make result variable
expInfo.cid = runCid;
expInfo.nneuron = nneuron;
expInfo.synIDs = synIDs;
expInfo.numExp = 2;
save(expInfoName,'expInfo');

%%Run experiments
for e = 1:length(activate)
    expFileName = sprintf('%sres%d.mat',cellDir,e);
    if 1;%~exist(expFileName,'file')

        nnRes.e = e;

        time1 = datetime('now');
        disp(sprintf('running experiment %d of %d',e,length(activate)))

        %%Define neuron experiment variables
        nneuron.pp{1}.AlphaSynapse = struct('node',activate{e},'gmax',0.01,'onset',20);% add an AlphaSynapse at this location
        nneuron.record{1}.cell.node = [1:numNodes];                            % record somatic voltage and voltage at synapse
        nneuron.record{1}.AlphaSynapse = struct('node',activate{e},'record','i');      % record synaptic current

        while 1
            pause(1)
            pass = 1;
            try
                out = t2n(nneuron,tree,'-w-q');              % execute t2n
            catch
                disp('neuron model failed. will try again')
                pass = 0;
            end
            if isempty(out.record{1})
                disp('record is empty. will try again')
                pass = 0;
            end
            if pass
                disp('model completed succesfully')
                break
            end
        end

        %% Collect Synapse Data
        synV = cat(2,out.record{1}.cell.v{uSynIDs})';
        nnRes.synV = synV(uSynIdx,:);

        %% Colect all node data
        allV = cat(2,out.record{1}.cell.v{:});
        maxV = max(allV,[],1);
        nnRes.maxV = maxV;

        if storeAllOut
            allOut(e).out = out;
        end

        %% Save results
        save(expFileName,'nnRes','-v7.3');

        %% Display time
        time2 = datetime('now');
        time3 = diff([time1 time2]);
        disp(sprintf('elapsed = %s',time3))


    end
end


%% Read in the result
load(expInfoName);
expInfo
numExp = expInfo.numExp;
%synVs = zeros(expNum,numSyn,1);
clear synVs
clear maxVs
%maxVs = zeros(expNum,numNodes);
c = 0;
for e = 1:numExp
    expFileName = sprintf('%sres%d.mat',cellDir,e);
    if exist(expFileName,'file')
        c = c+1;
        sprintf('found %s',expFileName)
        load(expFileName)
        synVs(c,:,:) =  nnRes.synV;
        maxVs(e,:) = nnRes.maxV;
    end
end
numRun = c;


%% plot the result (Vmem at soma and synapse and synaptic current)
clf
hold all
for s = 1:length(synIDs)
    plot(squeeze(synVs(1,s,:)));
end
legend('Synapse')
ylabel('Membrane potential [mV]')
% ylim([-85,0])
% xlim([0,50])
% subplot(2,1,2)
% plot(out.t,out.record{1}.AlphaSynapse.i{synIDs(end)})    % plot time vs synaptic current
% xlim([0,50])
% ylabel('Synaptic current [nA]')
% xlabel('Time [ms]')
pause(1)

%%
synVs(e,:,:) = synV(uSynIdx,:);
synMax = max(synVs,[],3);
synDev = synMax;
synDev = synDev - median(synDev(:));
synDev(synDev<0) = 0;
sumDev = sum(synDev,3);
allProp = sumDev;

for i = 1:numRun

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
pause(1)

%% Render with plot_tree

clf
props = maxVs;
props = props-min(props(:));
props = ceil(props * 255/(max(props(:))));
cmap = jet( 256 );
colormap(cmap)

for i = 1:numRun
    plot_tree(tree{1},props(i,:)');
    hold on
    scatter3(synPos(i,1),synPos(i,2),synPos(i,3),250,'k')
    %colorbar   % plot tree with region code
    set(gca,'clipping','off')
    pause(1)
end

%% polarity

influence = maxVs - neuron.params.v_init;
OFFbias = (influence(1,:)-influence(2,:))./(influence(1,:) + influence(2,:));

prop = OFFbias;
prop = prop + 1 * 100;
cmap = jet(200);
plot_tree(tree{1},prop');

%% Syn with Polarity
synVs(e,:,:) = synV(uSynIdx,:);
synMax = max(synVs,[],3);
synDev = synMax;
synDev = synDev - median(synDev(:));
synDev(synDev<0) = 0;
sumDev = sum(synDev,3);
allProp = sumDev;

meanV = squeeze(mean(mean(synVs,1),2));
difV = abs(meanV-median(meanV))>1;
showV = find(difV,1);
markerType = {'o','s'}
experimentFigureString = {'off', 'on'}

clf
hold on
axis 'equal'
axis 'off'
set(fig,'clipping','off')
hAx = gca;
set(hAx, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])

view(00,00)
s1 = scatter3(tr.X,tr.Y,tr.Z,1,'.','markeredgecolor',[.6 .6 .6]);
s3 = scatter3(tr.X(synIDs),tr.Y(synIDs),tr.Z(synIDs),100,'.','cdata',vColMap);

recordMovie = 1;
movieDir = sprintf('D:\\WorkDocs\\Publications\\VG3\\Revisions\\movies\\allVG3\\cid%d\\',runCid);
if ~exist(movieDir,'dir'),mkdir(movieDir);end
c = 0;

for r = 1


    for i = 1:numRun

        a = activate{i};
        eV = squeeze(synVs(i,:,:));
        sV = eV + 80;
        sV = ceil(sV*255/max(sV(:)));
        cmap = jet(256);
        tV = sV(:,1);
        vColMap = cmap(tV,:);
        s2 = scatter3(tr.X(a),tr.Y(a),tr.Z(a),100,markerType{i},'k','linewidth',1);
        title = experimentFigureString{i};

        for t = showV-10:showV+100
            t
            tV = sV(:,t);
            vColMap = cmap(tV,:);
            set(s3,'cdata',vColMap)
            % delete(s3)
            % s3 = scatter3(tr.X(synIDs),tr.Y(synIDs),tr.Z(synIDs),100,'.','cdata',vColMap);
            pause(.1)

            if recordMovie
                c = c+1;
                printName  = sprintf('%s%05.0f.png',movieDir,c);
                print(fig,printName,'-dpng')

            end
        end

        s2.delete;
     
    end

end

end




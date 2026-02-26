
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

%% Make default neuron

neuron = [];                                                % clear neuron structure
neuron.params.v_init = -60;                                 % starting membrane potential [mV] of all cells
neuron.params.dt = 0.025;                                   % integration time step [ms]
neuron.params.tstop = 300;                                  % stop simulation after this (simulation) time [ms]
neuron.params.prerun = -400;                                % add a pre runtime [ms] to let system settle
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
    tree{1}.R(1:2) = 3;     % make the first two nodes a new region
    tree{1}.rnames{3} = 'soma';   % name the region "soma"
    tree{1}.D(1:2) = 10;        % increase diameter of the somatic nodes
    tree = t2n_writeTrees(tree,[],mtrFile);                 % transform tree to NEURON morphology (.hoc file). This only has to be done once for each morphology
end
tree{1}.rnames{1} = 'dendrites';
tree{1}.rnames{2} = 'subtree';
tree{1}.Ri = 100;
tree{1}.Gm = 100 * 10^(-6);
tree{1}.Cm = 1;

%%Force fixed diameter
if 0
    tree{1}.D = tree{1}.D * 0 + median(tree{1}.D);
end

clf
plot_tree(tree{1},tree{1}.R);colorbar   % plot tree with region code
axis off

%% Get nep
nepFile = sprintf('%ssm2nrn_cid%d.mat',swcDir,runCid);
load(nepFile);
%load('Z:\Active\morganLab\MATLAB\TREES\morphos\sm2nrn_cid20.mat')

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


%% Make synapse experiment
nneuron = neuron;                                                         % copy standard neuron structure
plen = Pvec_tree(tree{1});                                                % get path length to soma at each node of morphology
%[~,synIDs] = max(plen);                                                   % search for the most far away point in the morpholgy
synIDs = nearest;


%%Force unique
[uSynIDs a uSynIdx] = unique(synIDs);


%%Make result variable
expInfo.cid = runCid;
expInfo.nneuron = nneuron;
expInfo.synIDs = synIDs;
expInfo.numExp = length(synIDs);
save(expInfoName,'expInfo');

for e = 1:length(synIDs)
    expFileName = sprintf('%sres%d.mat',cellDir,e);
    if ~exist(expFileName,'file')
        nnRes.e = e;

        time1 = datetime('now');
        disp(sprintf('running experiment %d of %d',e,numSyn))
        nneuron.pp{1}.AlphaSynapse = struct('node',synIDs(e),'gmax',0.01,'onset',20);% add an AlphaSynapse at this location
        if recordAllNodes
            nneuron.record{1}.cell.node = [1:numNodes];                            % record somatic voltage and voltage at synapse
        else
            nneuron.record{1}.cell.node = uSynIDs;                            % record somatic voltage and voltage at synapse
        end
        nneuron.record{1}.AlphaSynapse = struct('node',synIDs(e),'record','i');      % record synaptic current
        pause(1)

        while 1
            pass = 1

            try
                out = t2n(nneuron,tree,'-w-q');              % execute t2n
            catch
                disp('neuron model failed. will try again')
                pause(1)
                pass = 0
            end
            if isempty(out.record{1})
                disp('record is empty. will try again')
                pause(1)
                pass = 0;
            end
            if pass
                break
            end
        end

        synV = cat(2,out.record{1}.cell.v{uSynIDs})';
        if size(synV,2)>size(synVs,3), synVs(1,1,size(synV,2)) = 0;end
        %synVs(e,:,:) = synV(uSynIdx,:);
        nnRes.synV = synV(uSynIdx,:);


        if recordAllNodes
            allV = cat(2,out.record{1}.cell.v{:});
            maxV = max(allV,[],1);
            maxVs(e,:) = maxV;
            nnRes.maxV = maxV;

        end

        if storeAllOut
            allOut(e).out = out;
        end


        %% Save results
        save(expFileName,'nnRes','-v7.3');

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




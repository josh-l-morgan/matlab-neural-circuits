
%% Run NEURON from NEP with a range of conductances


%% SET Experiment parameters

testType = 3;
if testType == 1 % test conductances without inhibitory inputs
    runCids = [ 2 3 4 5 13 14 20]    
    testConductances =  10.^[-6:.5:-2];
    testInhibWeight = zeros(length(testConductances),1) + 10^4.5;
    runConThisTime = [1:length(testConductances)];
    showVideo = 0;
    expName = 'exOFFON_CondInhib_03';

elseif testType == 2
    
    runCids = [ 2 3 4 5 13 14 20]    
   
    eachCond = 10.^[-6:.5:-3];
    eachInhib = [10.^[-4:.5:-3]  10.^[-5 -4.5]]
    c = 0;
    clear testInhibWeight runConThisTime
    for tC = 1:length(eachCond)
        for tI = 1:length(eachInhib)
            c = c+1;
            testConductances(c) = eachCond(tC);
            testInhibWeight(c) = eachInhib(tI);
        end
    end
    runConThisTime = [1:c];
    showVideo = 0;
    expName = 'exOFFON_CondVsInhib_02';

elseif testType == 3 %test
       runCids = [ 4]    
   
    eachCond = 10.^[-6];
    eachInhib = 10^-4 .* [1/2 1/6 1/10];
  
    c = 0;
    clear testInhibWeight runConThisTime
    for tC = 1:length(eachCond)
        for tI = 1:length(eachInhib)
            c = c+1;
            testConductances(c) = eachCond(tC);
            testInhibWeight(c) = eachInhib(tI);
        end
    end
    runConThisTime = [1:c];
    expName = 'testOFFON_Cond08';
    showVideo = 0;
end
numCon = length(testConductances);


%% SET Cell parameters
restingPotential = -50;
cellReversalPotential = -50;
Ra = 100;
Cm = 1;
defaultConductance = 0.0003; % CHANGED IN EACH CONDITION

useExp2Syn = 1; %
withInhibit = 1;
useTimeCourse = 0;
onoffInhibScale = [1 .5];

%% SET Synapse parameters

startTime = 5;

syn1.e = 0; % Reversal in mV
syn2.e = -60; % Reversal in mV

%%Exp2Syn parameters

syn1.tau1 = 1;
syn1.tau2 = 1;%20;
syn1.i = 0;%defaultConductance;%26 * 10^-6;%0.01;
syn1.weight = 10^-4; % peak current in uS

syn2.tau1 = 1;
syn2.tau2 = 3;
syn2.i = 0;%defaultConductance;%9*10^-6;%0.01;
syn2.weight = 10^-2; % % peak current in uS, CHANGED IN EACH CONDITION

%%Alpha parameters
syn1.gmax = 10^-3;%26 * 10^-6; %0.001; % uS
syn1.onset = startTime;
%syn1.onsetSTD = 10;

syn2.gmax = 10^-4;
syn2.onset = startTime;
%syn2.onsetSTD = 1;

syn1.activeFrac = 1; %fraction of synapses activated each ms if timecourse is at max value
syn2.activeFrac = 1;

%% Set activation timecourse for clustering alphas

maxTime = 250;
etime = [0:1:maxTime];
rise1 = .5;
fall1 = .05;
rise2 = .05;
fall2 = 0.02;


y = 1-(1 * (1-rise1).^etime);
y2 =  (1 * (1-fall1).^etime)-1;
tc1 = y+y2;
tc1 = tc1/max(tc1);

y = 1-(1 * (1-rise2).^etime);
y2 =  (1 * (1-fall2).^etime)-1;
tc2 = y+y2;
tc2 = tc2/max(tc2);


cla, hold on
plot(tc1, 'g')
plot(tc2, 'r')
pause(.1)




%% Define neuron structure
neuron = [];                                                % clear neuron structure
neuron.params.v_init = restingPotential;                    % starting membrane potential [mV] of all cells
neuron.params.dt = 0.1;                                   % integration time step [ms] 0.025;
neuron.params.tstop = 15;                                  % stop simulation after this (simulation) time [ms]
neuron.params.prerun = -5; %400                           % add a pre runtime [ms] to let system settle
neuron.params.celsius = 33;                                 % temperature [celsius]
neuron.params.nseg = 'dlambda';                             % the dlambda rule is used to set the amount of segments per section. Alternatively, a fixed number can be entered or a string 'EachX' with X being the interval between segments [micron]
neuron.params.accuracy = 0;                                 % optional argument if number of segments should be increased in the soma and axon

%% Set neuron model mechanisms

neuron.experiment = 'test';                                                 % give your simulation run a name (necessary when using advanced t2n protocols)
for t = 1                                                                   % preconfiguration to loop through several morphologies, now only one morphology is considered
    neuron.mech{t}.all.pas = struct('g',defaultConductance,'Ra',Ra,'e',cellReversalPotential,'cm',Cm);    % add passive channel to all regions and define membrane capacity [µF/cm²], cytoplasmic resistivity [Ohm cm] and e_leak [mV]
    % neuron.mech{t}.soma.pas = struct('g',0.0004,'Ra',100,'e',-60,'cm',1);    % add passive channel to somatic regions (will overwrite the "all" definition) and define membrane capacity [µF/cm²], cytoplasmic resistivity [Ohm cm] and e_leak [mV]
    % neuron.mech{t}.all.k_ion.ek = -90;
    % neuron.mech{t}.all.na_ion.ena = 50;
    % neuron.mech{t}.soma.hh = struct('gnabar',0.25,'gkbar',0.036,'gl',0);     % add Hodgkin-Huxley Sodium channel only to soma
    neuron.record{t}.cell = struct('node',1,'record','v');                   % record voltage "v" from first node (i.e. the soma)
end


for v = 1:length(runCids)

    %% Define variables
    runCid = runCids(v);
    recordAllNodes = 1;
    storeAllOut = 0;
    w = warning ('on','all')

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
    ax = gca;

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
        tree = t2n_writeTrees(tree,[],mtrFile);                 % transform tree to NEURON morphology (.hoc file). This only has to be done once for each morphology
    end

    tree{1}.Ri = Ra;
    tree{1}.Gm = defaultConductance;
    tree{1}.Cm = Cm;

    %% Add artificial neurons
    trees = {};
    trees{1} = tree{1};
    trees{2} = struct('artificial','NetStim','start',startTime,'interval',15,'number',1); % add a NetStim as an artificial cell and define the start (10 ms) the interval (15 ms) and the number (10) of spikings
    trees{3} = struct('artificial','NetStim','start',startTime,'interval',15,'number',1); % add a NetStim as an artificial cell and define the start (10 ms) the interval (15 ms) and the number (10) of spikings
    %trees = t2n_writeTrees(trees,[],fullfile(pwd,'test.mtr'));                                    % tree morphologies are rewritten because this NetStim might be not written yet


    %% Fix diameters
    minDiam = 0.1;
    tree{1}.D(tree{1}.D<minDiam) = minDiam;
    setDiam = .5;
    if 0
        tree{1}.D = tree{1}.D * 0 + setDiam;
    end


    %% Show cell
    cla(ax)
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
    available = ones(numNodes,1);
    for s = 1:numSyn
        dists = sqrt((tr.X-synPos(s,1)).^2 + (tr.Y-synPos(s,2)).^2 + (tr.Z-synPos(s,3)).^2);
        dists(~available) = inf;
        minDists(s) = min(dists);
        close = find(dists == minDists(s),1);
        nearest(s) = close;
        available(close) = 0;
    end
    sharedNodes = length(nearest)-length(unique(nearest))

    cla(ax)
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
    amaCids = [];
    for i = 1:length(tis.cells.type.typeID)
        if tis.cells.type.typeID(i) == 7
            if intersect(useOnGroup,tis.cells.type.subTypeID(i))
                onBips = [onBips tis.cells.cids(i)];
            elseif intersect(useOffGroup,tis.cells.type.subTypeID(i))
                offBips = [offBips tis.cells.cids(i)];
            end
        elseif tis.cells.type.typeID(i) == 8
            if tis.cells.type.subTypeID(i) ~= 1; %not VG3
                amaCids = [amaCids tis.cells.cids(i)];
            end
        end
    end

    ONsyn = [];
    OFFsyn = [];
    AMCin = [];
    for i = 1:length(syn.pre);
        if intersect(onBips,syn.pre(i))
            ONsyn = [ONsyn i];
        elseif intersect(offBips,syn.pre(i))
            OFFsyn = [OFFsyn i];
        end
    end

    AMCin = find((syn.synType == 1) & ~(syn.pre==runCid));
    synIDs = nearest;
    activate{1} = synIDs(OFFsyn);
    activate{2} = synIDs(ONsyn);
    allInhibit = synIDs(AMCin);

    offBipPos = synPos(OFFsyn,:);
    onBipPos = synPos(ONsyn,:);


    clear inhibit
    inhibit{1} = [];
    inhibit{2} = [];

    if withInhibit
        if 0; % choose ON and OFF inhibitory synapses
            dThresh = 5;
            for d = 1:length(allInhibit)
                dists1 = sqrt((offBipPos(:,1)-synPos(AMCin(d),1)).^2 + (offBipPos(:,2)-synPos(AMCin(d),2)).^2 + ...
                    (offBipPos(:,3)-synPos(AMCin(d),3)).^2);
                dists2 = sqrt((onBipPos(:,1)-synPos(AMCin(d),1)).^2 + (onBipPos(:,2)-synPos(AMCin(d),2)).^2 + ...
                    (onBipPos(:,3)-synPos(AMCin(d),3)).^2);
                if sum(dists1<dThresh)
                    inhibit{1} = cat(1,inhibit{1},allInhibit(d));
                end
                if sum(dists2<dThresh)
                    inhibit{2} = cat(1,inhibit{2},allInhibit(d));
                end
            end

        else
            inhibit{1} = allInhibit;
            inhibit{2} = allInhibit;
        end
    end

    %% Show synapses
    cla(ax)
    hold on
    view(0,0)
    axis equal
    s1 = scatter3(tr.X,tr.Y,tr.Z,1,'.','markeredgecolor',[.6 .6 .6]);
    s2 = scatter3(tr.X(activate{1}),tr.Y(activate{1}),tr.Z(activate{1}),50,'o','markeredgecolor',[0 0 0.8],'LineWidth',2);
    s3 = scatter3(tr.X(activate{2}),tr.Y(activate{2}),tr.Z(activate{2}),50,'o','markeredgecolor',[1 0 0],'LineWidth',2);
    s4 = scatter3(tr.X(inhibit{1})+0.2,tr.Y(inhibit{1}),tr.Z(inhibit{1}),5,'o','markeredgecolor',[0 0 .8]);
    s5 = scatter3(tr.X(inhibit{2}),tr.Y(inhibit{2}),tr.Z(inhibit{2}),5,'o','markeredgecolor',[1 0 0]);
    set(gca,'clipping','off')


    %% Make synapse experiment

    plen = Pvec_tree(tree{1});                                                % get path length to soma at each node of morphology
    %[~,synIDs] = max(plen);                                                   % search for the most far away point in the morpholgy

    %%Force unique
    [uSynIDs a uSynIdx] = unique(synIDs);

    %%Make result variable
    expInfo.cid = runCid;
    expInfo.conductances = testConductances;
    expInfo.inhibWeight = testInhibWeight;
    expInfo.nneuron = neuron;
    expInfo.synIDs = synIDs;
    expInfo.numExp = 2;
    expInfo.expStr = {'activate OFF bipolar synapses','activate ON bipolar synapses'};
    expInfo.tree = tree{1};
    expInfo.OFFsyn = OFFsyn;
    expInfo.ONsyn = ONsyn;
    expInfo.AMCin = AMCin;

    save(expInfoName,'expInfo');

    for con = runConThisTime



        %% Run experiments
        for e = 1:length(activate)
            expFileName = sprintf('%sex%d_con%02.0f.mat',cellDir,e,con);
            %expFileName = sprintf('%sres%d.mat',cellDir,e);
            if 1%~exist(expFileName,'file')

                %%Refresh structures
                clear nnRes
                nnRes.e = e;
                nneuron = neuron;                                                         % copy standard neuron structure
                nnRes.nneuron = nneuron; %

                %% set test variables
                tree{1}.Gm = testConductances(con);
                nneuron.mech{1}.all.pas.g = testConductances(con);
                syn2.weight = testInhibWeight(con);

                time1 = datetime('now');
                condCount = find(runConThisTime==con);
                disp(sprintf('running cell %d of %d, condition %d of %d, experiment %d of %d',...
                    v,length(runCids),condCount,length(runConThisTime),e,length(activate)))



                %% Make synapses
                excite = activate{e};
                syn1.nodes = excite;
                syn2.nodes = inhibit{e};

                if useTimeCourse

                    rand1 = rand(length(syn1.nodes),length(tc1));
                    pass1 = rand1 <= repmat(tc1,[length(syn1.nodes) 1]);
                    rand2 = rand(length(syn1.nodes),length(tc1));
                    pass2 = pass1 & (rand2<= syn1.activeFrac);
                    [nID1 tID] = find(pass2>0);
                    onsets1 = etime(tID) + rand(size(tID));

                    rand1 = rand(length(syn2.nodes),length(tc2));
                    pass1 = rand1 <= repmat(tc2,[length(syn2.nodes) 1]);
                    rand2 = rand(length(syn2.nodes),length(tc2));
                    pass2 = pass1 & (rand2 <= syn2.activeFrac);
                    [nID2 tID] = find(pass2>0);
                    onsets2 = etime(tID) + rand(1,length(tID));
                else

                    nID1 = [1:length(syn1.nodes)];
                    onsets1 = nID1 * 0 + syn1.onset;
                    nID2 = [1:length(syn2.nodes)];
                    onsets2 = nID2 * 0 + syn2.onset;


                end


                %% Make Excitatory Synapses

                if useTimeCourse & ~useExp2Syn

                    exNodes = cat(1,syn1.nodes(nID1),syn2.nodes(nID2));
                    L1 = length(nID1);
                    L2 = length(nID2);
                    exE = cat(1,repmat(syn1.e,[L1 1]),repmat(syn2.e,[L2 1]));
                    exGmax = cat(1,repmat(syn1.gmax,[L1 1]),repmat(syn2.gmax * onoffInhibScale(e),[L2 1]));
                    exOnset = cat(1,repmat(syn1.onset,[L1 1]),repmat(syn2.onset,[L2 1]));

                    nneuron.pp{1}.AlphaSynapse = struct('node',exNodes,'gmax',exGmax,'onset',exOnset,'e',exE)
                    nneuron.record{1}.Exp2Syn = struct('node',exNodes,'record','i');           % record synaptic current

                else % no timecourse

                    exNodes = cat(1,syn1.nodes,syn2.nodes);
                    L1 = length(syn1.nodes);
                    L2 = length(syn2.nodes);
                    exE = cat(1,repmat(syn1.e,[L1 1]),repmat(syn2.e,[L2 1]));

                end
                if ~useExp2Syn %simple alpha synapses

                    exGmax = cat(1,repmat(syn1.gmax,[L1 1]),repmat(syn2.gmax * onoffInhibScale(e),[L2 1]));
                    exOnset = cat(1,repmat(syn1.onset,[L1 1]),repmat(syn2.onset,[L2 1]));
                    nneuron.pp{1}.AlphaSynapse = struct('node',exNodes,'gmax',exGmax,'onset',exOnset,'e',exE);% add an AlphaSynapse at this location
                    nneuron.record{1}.AlphaSynapse = struct('node',exNodes,'record','i');           % record synaptic current

                else
                    exTau1 = cat(1,repmat(syn1.tau1,[L1 1]),repmat(syn2.tau1,[L2 1]));
                    exTau2 = cat(1,repmat(syn1.tau2,[L1 1]),repmat(syn2.tau2,[L2 1]));
                    exI = cat(1,repmat(syn1.i,[L1 1]),repmat(syn2.i,[L2 1]));

                    exTags = {};
                    for g = 1:length(exNodes)
                        exTags{g} = ['t' num2str(exNodes(g))];
                    end
                    exTags = string(exTags);
                    exTag1 = exTags(1:L1);
                    exTag2 = exTags(L1+1:L1+L2);

                    nneuron.pp{1}.Exp2Syn =  struct('node',exNodes,'i',exI,'tau1',exTau1,'tau2',exTau2,'e',exE);% add an Exp2Syn at this location with 0.2 ms rise and 2.5 ms decay time
                    nneuron.record{1}.Exp2Syn = struct('node',exNodes,'record','i');           % record synaptic current

                    nneuron.con(1) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',syn1.nodes),'delay',0,'threshold',0.1,'weight',syn1.weight);  % connect the NetStim (cell 2) with the target (point process Exp2Syn of cell 1 at node specified in exIDs), and add threshold/weight and delay of the connection (NetStim parameters)
                    nneuron.con(2) = struct('source',struct('cell',3),'target',struct('cell',1,'pp','Exp2Syn','node',syn2.nodes),'delay',0,'threshold',0.1,'weight',syn2.weight * onoffInhibScale(e));  % connect the NetStim (cell 2) with the target (point process Exp2Syn of cell 1 at node specified in exIDs), and add threshold/weight and delay of the connection (NetStim parameters)


                end


                %% Record all nodes
                nneuron.record{1}.cell.node = [1:numNodes];                            % record somatic voltage and voltage at synapse
                % nneuron.pp{1}.SEClamp = struct('node',1,'times',[10],'amp',[0]);   
                % nneuron.record{1}.SEClamp = struct('node',1,'record','i');                      % record current "i" from SEClamp electrode which will be located at node 1



                %% Run Model
                while 1
                    pause(1)
                    pass = 1;
                    try
                        out = t2n(nneuron,trees,'-w-q');            % execute t2n
                    catch ME
                        disp(['Error = ' ME.message])
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
                pause(.01)
              
                figure(fig)

                %% Collect Synapse Data
                synV = cat(2,out.record{1}.cell.v{uSynIDs})';
                nnRes.synV = synV(uSynIdx,:);

                try synI = cat(2,out.record{1}.Exp2Syn.i{exNodes});

                catch
                    try synI = cat(2,out.record{1}.AlphaSynapse.i{syn1.nodes});
                    catch
                        synI = [];
                    end
                end
                nnRes.synI = synI;

                %% Colect all node data
                allV = cat(2,out.record{1}.cell.v{:});
                maxV = max(allV,[],1);
                nnRes.maxV = maxV;
                %nnRes.allV = allV;

                if storeAllOut
                    allOut(e).out = out;
                end

                %% Save results

                nnRes.time = out.t;
                nnRes.rootV = out.record{1}.cell.v{1};
                nnRes.syn1 = syn1;
                nnRes.syn2 = syn2;
                save(expFileName,'nnRes','-v7.3');

                %% Display time
                time2 = datetime('now');
                time3 = diff([time1 time2]);
                disp(sprintf('elapsed = %s',time3))


            end
        end

        set(fig,'visible','on')

        %% Read in the result
        load(expInfoName);
        numExp = expInfo.numExp;
        %synVs = zeros(expNum,numSyn,1);
        clear synVs maxVs synIs rootV
       
        %allVs = zeros(numExp,size(nnRes.allV,2),size(nnRes.allV,1));
        %maxVs = zeros(expNum,numNodes);
        c = 0;
        for e = 1:numExp
            expFileName = sprintf('%sex%d_con%02.0f.mat',cellDir,e,con);
            if exist(expFileName,'file')
                c = c+1;
                % disp(sprintf('found %s',expFileName))
                load(expFileName);
                synVs(c,:,:) =  nnRes.synV;
                synIs{c} = nnRes.synI;
                maxVs(e,:) = nnRes.maxV;
                rootV(e,:) = nnRes.rootV;
                %allVs(c,:,:) = nnRes.allV';
            end
        end
        numRun = c;
        nTime = nnRes.time;


        %% plot the result (Vmem at soma and synapse and synaptic current)

        figure(f)
        clf
        set(fig,'Visible','on');
        

        sp = subplot(4,2,1); hold on
        plot(sp,nTime,rootV(1,:),'b')
        plot(sp,nTime,rootV(2,:),'r')
        xlim([0 max(nTime)])

        sp = subplot(4,2,3); hold on
        
        synI = synIs{1};
        for n = 1:size(synI,2)
            plot(sp,nTime,synI(:,n)','b','linewidth',.01)
        end

        synI = synIs{2};
        for n = 1:size(synI,2)
            plot(sp,nTime+1,synI(:,n)','r','linewidth',.01)
        end
        ylabel({'current'})
        xlim([0 max(nTime)])


        sp = subplot(4,2,2); cla(sp)
        hold on
        nodes = nnRes.syn1.nodes;
        for n = 1:length(expInfo.OFFsyn);
            plot(sp,nTime,squeeze(synVs(1,expInfo.OFFsyn(n),:)),'b','linewidth',.01)
        end
        %plot(nTime,meanV,'b','lineWidth',2)


        nodes = nnRes.syn1.nodes;
        for n = 1:length(expInfo.OFFsyn);
            plot(sp,nTime+1,squeeze(synVs(2,expInfo.OFFsyn(n),:)),'r','linewidth',.01)
        end
        ylabel({'offSyn'})
        xlim([0 max(nTime)])


        sp = subplot(4,2,4); cla(sp),
        hold on

        nodes = nnRes.syn2.nodes;
        for n = 1:length(expInfo.ONsyn);
            plot(sp,nTime,squeeze(synVs(1,expInfo.ONsyn(n),:)),'b','linewidth',.01)
        end
        %plot(nTime,meanV,'r','lineWidth',2)
        %ylim([-70 -20])
        %xlim([8 15])

        nodes = nnRes.syn2.nodes;
        for n = 1:length(expInfo.ONsyn);
            plot(sp,nTime+1,squeeze(synVs(2,expInfo.ONsyn(n),:)),'r','linewidth',.01)
        end
        %ylim([-70 20])
        %xlim([8 150])

        ylabel({'onSyn'})
        xlim([0 max(nTime)])

        pause(.1)

        if 0 %if all voltages are recored
            cla
            subplot(2,1,1,ax),   hold on
            nodes = nnRes.syn1.nodes;
            meanV = squeeze(mean(allVs(1,nodes,:),2));
            for n = 1:length(nodes);
                plot(ax,nTime,squeeze(allVs(1,nodes(n),:)),'b','linewidth',.01)
            end
            %plot(nTime,meanV,'b','lineWidth',2)

            nodes = nnRes.syn2.nodes;
            meanV = squeeze(mean(allVs(1,nodes,:),2));
            for n = 1:length(nodes);
                plot(ax,nTime+1,squeeze(allVs(1,nodes(n),:)),'r','linewidth',.01)
            end
            %plot(nTime,meanV,'r','lineWidth',2)
            ylim([-70 -20])
            %xlim([8 15])


            subplot(2,1,2,ax),  hold on
            nodes = nnRes.syn1.nodes;
            for n = 1:length(nodes);
                plot(ax,nTime,squeeze(allVs(2,nodes(n),:)),'b','linewidth',.01)
            end

            nodes = nnRes.syn2.nodes;
            for n = 1:length(nodes);
                plot(ax,nTime+1,squeeze(allVs(2,nodes(n),:)),'r','linewidth',.01)
            end
            ylim(ax,[-70 20])
            %xlim([8 150])

            pause(.1)

        end


        %% color synapses
        if 0
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
                cla
                scatter3(tr.X,tr.Y,tr.Z,'.','k')
                axis 'equal'
                hold on
                synScat = scatter3(synPos(:,1),synPos(:,2),synPos(:,3));
                set(synScat,'CData',prop)
                scatter3(tr.X(nearest),tr.Y(nearest),tr.Z(nearest),'p','g')

            end
            pause(.1)
        end

        %% Render voltage with plot_tree
        if 0
            cla
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
                pause(.1)
            end

        end

        %% polarity
         figure(f)
        set(fig,'Visible','on')
        sp = subplot(2,1,2);
        cla(sp)
        cmap = jet(200);
        colormap(cmap)

        OFFbias = (maxVs(2,:)-maxVs(1,:))./sum(maxVs,1);
        medBias = median(abs(OFFbias));
        prop = OFFbias;
        span = 50/medBias;
        prop = prop * span + 100;
        prop = round(prop);
        prop(prop<1) = 1;
        prop(prop>200) = 200;
        if sum(isnan(prop))
            disp('WARNING: Some results are not numbers');
        end
        prop(isnan(prop)) = 100;
        sc = scatter3(sp,tr.X,tr.Y,tr.Z,5,'o','k','markerfacecolor','flat');
        set(sc,'cdata',cmap(prop,:))
        axis equal
        view(2,2)

        drawnow


        %% Syn with Polarity
        if showVideo

            synVs(e,:,:) = synV(uSynIdx,:);
            synMax = max(synVs,[],3);
            synDev = synMax;
            synDev = synDev - median(synDev(:));
            synDev(synDev<0) = 0;
            sumDev = sum(synDev,3);
            allProp = sumDev;

            meanV = squeeze(mean(mean(synVs,1),2));
            difV = abs(meanV-median(meanV))>.01;
            showV = find(difV,1);
            markerType = {'o','o','d'}
            experimentFigureString = {'off', 'on'}

            cla
            hold on
            axis 'equal'
            axis 'off'
            set(fig,'clipping','off')
            hAx = gca;
            set(hAx, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])

            view(00,00)
            s1 = scatter3(tr.X,tr.Y,tr.Z,1,'.','markeredgecolor',[.6 .6 .6]);
            % s1 = scatter3(tr.X,tr.Y,tr.Z,10,'o','markerfacecolor','flat','markerfacealpha',.2,...
            %     'markeredgecolor','none');
            s4 = scatter3(tr.X(nnRes.syn2.nodes),tr.Y(nnRes.syn2.nodes),tr.Z(nnRes.syn2.nodes),...
                50,markerType{3},'k','linewidth',2);

            s3 = scatter3(tr.X(synIDs),tr.Y(synIDs),tr.Z(synIDs),50,'o','markerfacecolor','flat','markerfacealpha',1,...
                'markeredgecolor','none');

            recordMovie = 0;
            movieDir = sprintf('D:\\WorkDocs\\Publications\\VG3\\Revisions\\movies\\allVG3\\cid%d\\',runCid);
            if ~exist(movieDir,'dir'),mkdir(movieDir);end
            c = 0;

            for r = 1


                for i = 1:numRun

                    a = activate{i};
                    eV = squeeze(synVs(i,:,:));
                    eSum = sum(abs(eV-eV(1,1)),1);
                    eStart = find(eSum>0,1)
                    %eV = squeeze(allVs(i,:,:));
                    sV = eV + 90;
                    sV = ceil(sV*255/max(sV(:)));
                    sV(sV<1)=1;
                    cmap = jet(256);
                    tV = sV(:,1);
                    vColMap = cmap(tV,:);
                    s2 = scatter3(tr.X(a),tr.Y(a),tr.Z(a),100,markerType{i},'k','linewidth',2);
                    title = experimentFigureString{i};


                    for t = eStart:eStart+200;%size(sV,2)
                        t
                        tV = sV(:,t);
                        vColMap = cmap(tV,:);
                        set(s3,'cdata',vColMap)
                        %set(s1,'cdata',vColMap)

                        % delete(s3)
                        % s3 = scatter3(tr.X(synIDs),tr.Y(synIDs),tr.Z(synIDs),100,'.','cdata',vColMap);
                        pause(.01)

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

    end
end

disp('finished running models')

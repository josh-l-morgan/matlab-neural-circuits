function[s] = parseExcelSynProps(dat)
%%Dat should be copy and paste excel sheet with syn data

%% Clear dat
if 0
    dat = {};
    %%Paste KxR cell Registry -> SynapseTyping
end
[datH datW] = size(dat);

s.emRes = [.004 .004 .040];

%% parse headers
headers = dat(1,:);
colNames = {'cid #', 'cidCol';'branch point', 'bpCol'; 'clearRGC', 'clearRGCCol';...
    'large','largeCol'; 'glia','gliaCol'; 'spine', 'spineCol';...
    'large vec', 'vecCol'; 'light mito','mitoCol'; 'input','inputCol';...
    'pre cid','preCidCol'; 'synapse point', 'synPointCol'; 'syn type','synTypeCol';...
    'junction point', 'jpCol'; 'spine point', 'spinePointCol'};

%%Create empty column variables
for c = 1:length(colNames)
    eval(sprintf('%s = [];',colNames{c,2}));
end
%%Find columns
for i = 1:datW
    head = dat{1,i};
    for c = 1:length(colNames);
        if strcmp(head,colNames{c,1})
            eval(sprintf('%s = %d;',colNames{c,2},i));
        end
    end
end

%% Parse post cid
currCid = 0;
s.postCid = zeros(datH,1);
for i = 1:datH
    val = dat{i,cidCol};
    if ~strcmp(class(val),'char')
        if ~isempty(val)
            currCid = val;
        end
    end
    s.postCid(i) = currCid;
end

%% Parse pre cid
s.preCid = zeros(datH,1);
for i = 1:datH
    val = dat{i,cidCol};
    if ~strcmp(class(val),'char')
        if ~isempty(val)
            s.preCid(i) = val;
        end
    end
end

%% Parse syn Type
s.synTypeNames = {'rgc','unk','lin','ctx','trn','not','?'};
s.synType = zeros(datH,1);
for h = 1:datH
    val = dat{h,synTypeCol};
    if strcmp(class(val),'char')
        for c = 1:length(synTypeNames)
            if sum(regexp(val,synTypeNames{c}))
                s.synType(h,1) = c;
                break
            end
        end
    end
end

%% Parse syn Prop
s.hasSynProp = zeros(datH,1)-1;
s.clearRGC = zeros(datH,1)-1;
s.large = zeros(datH,1)-1;
s.hasGlia = zeros(datH,1)-1;
s.hasSpine = zeros(datH,1)-1;
s.largeVec = zeros(datH,1)-1;
s.lightMito = zeros(datH,1)-1;
s.hasInput = zeros(datH,1)-1;

for h = 2:datH
    if ~isempty(dat{h,clearRGCCol})
        s.clearRGC(h,1) = double(dat{h,clearRGCCol})';
    end
     if ~isempty(dat{h,largeCol})
        s.large(h,1) = double([dat{h,largeCol}])';
    end
    if ~isempty(dat{h,gliaCol})
        s.hasGlia(h,1) = double([dat{h,gliaCol}])';
        s.hasSynProp(h,1) = 1;
    end
    if ~isempty(dat{h,spineCol})
        s.hasSpine(h,1) = double([dat{h,spineCol}])';
    end
    if ~isempty(dat{h,vecCol})
        s.largeVec(h,1) = double([dat{h,vecCol}])';
    end
    if~isempty(dat{h,mitoCol})
        s.lightMito(h,1) = double([dat{h,mitoCol}])';
    end
    if ~isempty(dat{h,inputCol})
        s.hasInput(h,1) = double([dat{h,inputCol}])';
    end
end

%% Parse Pos
s.isPoint = zeros(datH,1); %row is a point recorded
s.pos = zeros(datH,3);
for h = 1:datH
    val = dat{h,synPointCol};
    if ~isempty(val)
        anch = sscanf(val,'(%d, %d, %d');
        if size(anch,1) == 3
            s.isPoint(h,1) = 1;
            s.pos(h,1:3) = [anch(1) * emRes(1) anch(2) * emRes(2) anch(3) * emRes(3)];
        end
    end
end






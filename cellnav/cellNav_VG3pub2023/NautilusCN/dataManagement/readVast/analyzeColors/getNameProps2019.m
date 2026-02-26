function[nameProps] = getNameProps(objNames,aliases)


nameProps.oldNames = objNames;


%% apply aliases
load('MPN.mat')
if ~exist('aliases','var')
    aliases = [];
    if exist([MPN 'dat.mat'],'file')
        load([MPN 'dat.mat']);
    end
    try
        aliases = dat.alias;
    end
end


nams = objNames;
for i = 1:length(objNames)
    nams{i}  = namAliases(objNames{i},aliases);
end
nameProps.names = nams;


%% target string

if ~exist('targString','var')
    targString = 'no target string was entered by user';
end
targString = lower(targString);

%% parse cids
allIds = cell(length(nams));
for i = 1:length(nams)
    
    nam = lower(nams{i});
    nameProps.parseError(i) = 0;
    nameProps.doubt(i) = sum(regexp(nam,'?')>0);
    
    %%Get IDs
    [cids cidPos parseError] = getCids(nam);
    nameProps.cids{i} = cids;
    nameProps.cidPos{i} = cidPos;
    nameProps.parseError(i) = parseError;
    
    %%legacy ids
    allCellIds = cids;
    tag.allIDs{i} = allCellIds;
    if isempty(allCellIds)
        nameProps.cellNum(i) = 0; % no cell number found
    else
        nameProps.cellNum(i) = allCellIds(1);
    end
    
end



%% parse tags
for i = 1:length(nams)
    nam = lower(nams{i});
    nam = namAliases(nam,aliases);
    tag.firstNum(i) = getFirstNumber(nam);
    tag.cell(i) = sum(regexp(nam,'cell')>0);
    tag.hasCid(i) = sum(regexp(nam,'cid')>0);
    tag.part(i) = sum(regexp(nam,'part')>0);
    tag.frag(i) = sum(regexp(nam,'frag')>0);
    tag.indifrag(i) = sum(regexp(nam,'indifrag')>0);
    tag.seed(i) = sum(regexp(nam,'seed')>0);
    tag.somatic(i) = sum(regexp(nam,'somatic')>0);
    tag.soma(i) = sum((regexp(nam,'soma')>0)) + sum((regexp(nam,'cb')>0));
    tag.sft(i) = sum(regexp(nam,'sft')>0);
    
    tag.vgc(i) = sum(regexp(nam,'vgc')>0);
    tag.bpc(i) = sum(regexp(nam,'bpc')>0);
    tag.rgc(i) = sum(regexp(nam,'rgc')>0);
    tag.amc(i) = sum(regexp(nam,'amc')>0);
    tag.glia(i) = sum(regexp(nam,'glia')>0);
    tag.mgc(i) = sum(regexp(nam,'mgc')>0);
    
    tag.ribbon(i) = sum(regexp(nam,'ribbon')>0);
    tag.vec(i) = sum(regexp(nam,'vec')>0);
    tag.mito(i) = sum(regexp(nam,'mito')>0);
    tag.cent(i) = sum(regexp(nam,'cent')>0);
    tag.nuc(i) = sum(regexp(nam,'nuc')>0);
    tag.nucl(i) = sum(regexp(nam,'nucl')>0);
    tag.golg(i) = sum(regexp(nam,'golg')>0);
    tag.hchr(i) = sum(regexp(nam,'hchr  ')>0);
    tag.er(i) = sum(regexp(nam,'er')>0);
    tag.rer(i) = sum(regexp(nam,'rer')>0);
    tag.prc(i) = sum(regexp(nam,'prc')>0);
    tag.rbs(i) = sum(regexp(nam,'rbs')>0);
    tag.vac(i) = sum(regexp(nam,'vac')>0);
    tag.drkstk(i) = sum(regexp(nam,'drkstk')>0);
    tag.brh(i) = sum(regexp(nam,'brh')>0);
    tag.init(i) = sum(regexp(nam,'init')>0);
    
    tag.syn(i) = sum(regexp(nam,'syn')>0);
    tag.dyad(i) = sum(regexp(nam,'dyad')>0);
    tag.rib(i) = sum(regexp(nam,'rib')>0);
    tag.gap(i) = sum(regexp(nam,'gap')>0);
    tag.jnc(i) = sum(regexp(nam,'jnc')>0);
    
    
    
    tag.LMB(i) = sum(regexp(nam,'lmb')>0);
    tag.DMB(i)  = sum(regexp(nam,'dmb')>0);
    tag.ldm(i) = sum(regexp(nam,'ldm')>0);
    tag.sdm(i) = sum(regexp(nam,'sdm')>0);
    tag.flj(i) = sum(regexp(nam,'flj')>0);
    tag.junc(i) = sum(regexp(nam,'junc')>0);
    tag.tcr(i) = sum(regexp(nam,'tcr')>0);
    tag.rgc(i) = sum(regexp(nam,'rgc')>0);
    tag.lin(i) = sum(regexp(nam,' lin')>0);
    tag.axprim(i) = sum(regexp(nam,'axprim')>0);
    tag.axon(i) = sum(regexp(nam,'axon')>0);
    tag.unk(i) = sum(regexp(nam,'unk')>0);
    tag.isp(i) = sum(regexp(nam,'isp')>0);
    tag.distal(i) = sum(regexp(nam,'distal')>0);
    
    tag.targS(i) = sum(regexp(nam,targString)>0);
    tag.bouton(i) = sum(regexp(nam,'bouton')>0);
    tag.bout(i) = sum(regexp(nam,'bout ')>0);
    tag.glom(i) = sum(regexp(nam,'glom')>0);
    tag.glm(i) = sum(regexp(nam,'glm')>0);
    tag.center(i) = sum(regexp(nam,'center')>0);
    tag.retZone(i) = sum(regexp(nam,'retzone')>0);
    tag.toamc(i) = sum(regexp(nam,'tamc')>0); %Karl 05/23/22 - changed toamc to tamc
    tag.torgc(i) = sum(regexp(nam,'trgc')>0); %Karl 05/23/22 - changed torgc to trgc
    
end

nameProps.ofID = tag.indifrag | tag.init |tag.brh |tag.frag | tag.part ...
    | tag.cell | tag.glom | tag.center | tag.axprim |tag.retZone |...
    tag.distal | tag.bout | tag.soma | tag.seed;
nameProps.tag = tag;


%% Parse edges

syns = find(tag.syn & ~nameProps.doubt);
ribs = find(tag.rib & ~nameProps.doubt);
dyads = find(tag.dyad & ~nameProps.doubt);
torgc = find(tag.torgc & ~nameProps.doubt); %Karl 05/23/22 - added lines
toamc = find(tag.toamc & ~nameProps.doubt); %Karl 05/23/22 - added lines

edgeSegs = cat(2,syns,ribs,dyads);
edgeTypes = cat(2,syns*0+1,ribs*0+2,dyads*0+3,torgc * 0 +4,toamc * 0 + 5);
edgeTags = {' syn ' ' rib ' ' dyad '};

edges = [];
postGroups = [];
nameProps.synProp = {};
for i = 1:length(edgeSegs)
    nam = nams{edgeSegs(i)};
    [conProp] = parseEdgeSurround(nam,edgeTags{edgeTypes(i)});
    if ~isempty(conProp.edges)
        conProp.segID = edgeSegs(i);
        edge1 = conProp.edges(:,[2 1]);
        edge1(:,3) = edgeSegs(i);
        edge1(:,4) = edgeTypes(i);
        edges = cat(1,edges,edge1);
        postGroups = cat(1,postGroups,conProp.postGroup);
    end
    nameProps.synProp{i} = conProp;
    
end
nameProps.edges = edges;

%% parse junctions
jncs = find((tag.jnc ) & ~nameProps.doubt);
juncs = find((tag.junc) & ~nameProps.doubt);

edgeSegs = cat(2,jncs,juncs);
edgeTypes = cat(2,jncs*0+1,juncs*0+2);
edgeTags = {'jnc' 'junc'};

jncEdge = [];
for i = 1:length(edgeSegs)
    nam = nams{edgeSegs(i)};
    [conProp] = parseEdgeSurround(nam,edgeTags{edgeTypes(i)});
    edge1(:,3) = edgeSegs(i);
    jncEdge = cat(1,jncEdge,edge1);
    nameProps.juncProps{i} = conProp;
end

nameProps.juncs = jncEdge;



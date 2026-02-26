function[nameProps] = getNameProps(objNames)


nameProps.names = objNames;

if ~exist('targString','var')
    targString = 'no target string was entered by user';
end


load('MPN.mat')
if exist([MPN 'dat.mat'],'file')
    load([MPN 'dat.mat']);
    try, aliases = dat.alias; end
end
    


targString = lower(targString);

allIds = cell(length(objNames));
for i = 1:length(objNames)    
    
    nam = lower(objNames{i});
    nameProps.parseError(i) = 0;
    nameProps.doubt(i) = sum(regexp(nam,'?')>0);

    %%Get IDs
    [cids cidPos parseError] = getCids(nam);
    if exist('aliases','var')
        cids = checkAlias(cids,aliases);
    end
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




for i = 1:length(objNames)  
    nam = lower(objNames{i});
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
       
end

nameProps.ofID = tag.indifrag | tag.init |tag.brh |tag.frag | tag.part ...
    | tag.cell | tag.glom | tag.center | tag.axprim |tag.retZone |...
    tag.distal | tag.bout | tag.soma | tag.seed;
nameProps.tag = tag;


%% Parse edges

syns = find(tag.syn & ~nameProps.doubt);
ribs = find(tag.rib & ~nameProps.doubt);
dyads = find(tag.dyad & ~nameProps.doubt);

edgeSegs = cat(2,syns,ribs,dyads);
edgeTypes = cat(2,syns*0+1,ribs*0+2,dyads*0+3);
edgeTags = {' syn ' ' rib ' ' dyad '};

edges = [];
postGroups = [];
nameProps.synProp = {};
for i = 1:length(edgeSegs)
    nam = objNames{edgeSegs(i)};
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
    nam = objNames{edgeSegs(i)};
    [conProp] = parseEdgeSurround(nam,edgeTags{edgeTypes(i)});
    edge1(:,3) = edgeSegs(i);
    jncEdge = cat(1,jncEdge,edge1);
    nameProps.juncProps{i} = conProp; 
end

nameProps.juncs = jncEdge;



%% Legacy synapse
%{
%% Get edges
%edges
numEdge = 0;
numJuncs = 0;
numGlm = 0;
edges = [];
juncs = [];
synType = [];
synProp = [];
 numIsp = 0;
    inSpine = [];
for i = 1:length(nameProps.names)
    if nameProps.part(i)
       parts = nameProps.allIDs{i};
       if length(parts) >1
                numEdge = numEdge+1;
  
       edges(numEdge,:) = [parts(1:2) i];
       end
    end
    
    if nameProps.glm(i)
        parts = nameProps.allIDs{i};
        if length(parts) >1
            numGlm = numGlm+1;
            nam = nameProps.names{i};
            dlim = [0 regexp(nam,' ') length(nam)+1];
            clear boutonNum spineNum
            for d = 1:length(dlim)-1
                phrase = nam(dlim(d)+1:dlim(d+1)-1);
                gPos = regexp(phrase,'glm');
               
                if sum(gPos);
                    boutonNum = str2num(phrase(1:gPos-1));
                    suffix = phrase(gPos+3:end);
                    filimentous = double(sum(regexp(suffix,'F'))>0);
                    axoAxonic = double(sum(regexp(suffix,'B'))>0);
                    inhibitory = double(sum(regexp(suffix,'I'))>0);
                    isDend = sum(regexp(suffix,'D'))>0;
                    
                    spinePos = regexp(suffix,'S');
                    if ~sum(spinePos)
                        spineNum = 0;
                    else
                        afterS = [];
                        for sp = spinePos+1:length(suffix);
                            
                            if isempty(str2num(suffix(sp)))
                                break
                            else
                                afterS = [afterS suffix(sp)];
                            end
                        end
                        
                        if isempty(afterS)
                            spineNum = 1;
                        else
                            spineNum = str2num(afterS);
                        end
                    end
                
                synProp(numGlm).pre = parts(2);
                synProp(numGlm).post = parts(1);
                synProp(numGlm).objectID = i;
                synProp(numGlm).isDend = isDend;
                synProp(numGlm).boutonNeighbor = boutonNum;
                synProp(numGlm).spineNum = spineNum;
                synProp(numGlm).filimentous = filimentous;
                synProp(numGlm).axoAxonic = axoAxonic;
                synProp(numGlm).lin = axoAxonic | inhibitory;
                synProp(numGlm).tcr = ~(axoAxonic | inhibitory);
                break
                end
                
                
            end
        end
    end
    
   
   if nameProps.isp(i)
        parts = nameProps.allIDs{i};
        if length(parts) >1
            numGlm = numGlm+1;
            nam = nameProps.names{i};
            dlim = [0 regexp(nam,' ') length(nam)+1];
            clear boutonNum spineNum
            for d = 1:length(dlim)-1
                phrase = nam(dlim(d)+1:dlim(d+1)-1);
                gPos = regexp(phrase,'isp');
               
                if sum(gPos);
                    numIsp = numIsp+1

                    inNum = str2num(phrase(1:gPos-1));
                    tipNum = str2num(phrase(gPos+3:end));
                    if isempty(inNum)
                        inNum = 0;
                    end
                    if isempty(tipNum)
                        tipNum = 0;
                    end
                
                inSpine(numIsp).in = inNum;
                inSpine(numIsp).tip = inNum;
                inSpine(numIsp).pre = parts(2);
                inSpine(numIsp).post = parts(1);
                inSpine(numIsp).objectID = i;
                
                break
                end
                
                
            end
        end
    end
           
    
     if nameProps.junc(i)
       parts = nameProps.allIDs{i};
       if length(parts) >1
               numJuncs = numJuncs+1;
   
       juncs(numJuncs,:) = [parts(1:2) i];
       end
     end
end
% 
% postList = sort(unique(edges(:,1)));
% preList = sort(unique(edges(:,2)));

nameProps.edges = edges;
nameProps.juncs = juncs;
nameProps.synProp = synProp;
nameProps.inSpine = inSpine;

%% Make connectivity matrix.
if isempty(edges)
    preID = [];
    postID = [];
else
preID = unique(edges(:,1));
postID = unique(edges(:,2));
end

%}




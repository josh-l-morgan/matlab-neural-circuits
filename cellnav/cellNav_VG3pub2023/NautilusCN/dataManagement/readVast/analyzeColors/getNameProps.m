function[nameProps] = getNameProps(objNames,targString)

if ~exist('targString','var')
    targString = 'no target string was entered by user';
end

targString = lower(targString);

allIds = cell(length(objNames));
for i = 1:length(objNames)    
    postNum{i} = ' ';
    nam = lower(objNames{i});
    ascNam = uint8(nam);
    isNum = (ascNam>=48) & (ascNam<=57);
    
    %%Get cell id numbers
    tempNum = [];
    allCellIds = [];
    for s = 1:length(nam)
        
            if isNum(s)
                tempNum = [tempNum nam(s)];
            else
                if ~isempty(tempNum)
                    allCellIds = [allCellIds str2num(tempNum)];
                    postNum{i} = nam(s); 
                    tempNum = [];
                end
            end
    end
    if ~isempty(tempNum)
         allCellIds = [allCellIds str2num(tempNum)];
        postNum{i} = nam(s); 
          tempNum = [];
    end
    
    allIDs{i} = allCellIds;
    if isempty(allCellIds)
        cellNum(i) = 0; % no cell number found
    else
        cellNum(i) = allCellIds(1);
    end
    
    isCell(i) = sum(regexp(nam,'cell'));
    part(i) = sum(regexp(nam,'part'));
    LMB(i) = sum(regexp(nam,'lmb'));
    DMB(i)  = sum(regexp(nam,'dmb'));
    ldm(i) = sum(regexp(nam,'ldm'));
    sdm(i) = sum(regexp(nam,'sdm'));
    syn(i) = sum(regexp(nam,'syn'));
    flj(i) = sum(regexp(nam,'flj'));
    junc(i) = sum(regexp(nam,'junc'));
    frag(i) = sum(regexp(nam,'frag'));
    tcr(i) = sum(regexp(nam,'tcr'));
    rgc(i) = sum(regexp(nam,'rgc'));
    lin(i) = sum(regexp(nam,' lin'));
    axprim(i) = sum(regexp(nam,'axprim'));
    axon(i) = sum(regexp(nam,'axon'));
    unk(i) = sum(regexp(nam,'unk'));
    isp(i) = sum(regexp(nam,'isp'));
    distal(i) = sum(regexp(nam,'distal'));
   
    targS(i) = sum(regexp(nam,targString));
    bouton(i) = sum(regexp(nam,'bouton'));
    bout(i) = sum(regexp(nam,'bout '));
    glom(i) = sum(regexp(nam,'glom'));
    glm(i) = sum(regexp(nam,'glm'));
    center(i) = sum(regexp(nam,'center'));
    retZone(i) = sum(regexp(nam,'retzone'));
    
    
end

ofID = frag | part | isCell | glom | center | axprim |retZone | distal | bout;

%% place in structure
nameProps.names = objNames;
nameProps.cellNum = cellNum;
nameProps.allIDs = allIDs;
nameProps.tcr = tcr>0;
nameProps.lin = lin>0;
nameProps.rgc = rgc>0;
nameProps.axon = axon>0;
nameProps.retZone = retZone;
nameProps.postNum = postNum;
nameProps.cell = isCell;
nameProps.part = part;
nameProps.ldm = (ldm>0) & (tcr==0);
nameProps.sdm = sdm>0;
nameProps.unk = unk>0;
nameProps.LMB = LMB;
nameProps.DMB = DMB;
nameProps.syn = syn;
nameProps.flj = flj;
nameProps.junc = junc;
nameProps.targS = targS;
nameProps.bouton = bouton;
nameProps.bout = bout;
nameProps.frag = frag;
nameProps.ofID = ofID;
nameProps.glom = glom;
nameProps.glm = glm;
nameProps.isp = isp;


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




function[nameProps] = getMicroNameProps(objNames,targString)

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
    glom(i) = sum(regexp(nam,'glom'));
    glm(i) = sum(regexp(nam,'glm'));
    center(i) = sum(regexp(nam,'center'));
    retZone(i) = sum(regexp(nam,'retzone'));
    
    
end

ofID = frag | part | isCell | glom | center | axprim |retZone | distal;

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





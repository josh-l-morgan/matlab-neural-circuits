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
    syn(i) = sum(regexp(nam,'syn'));
    flj(i) = sum(regexp(nam,'flj'));
    junc(i) = sum(regexp(nam,'junc'));
    frag(i) = sum(regexp(nam,'frag'));
    tcr(i) = sum(regexp(nam,'tcr'));
    lin(i) = sum(regexp(nam,' lin'));
   
   targS(i) = sum(regexp(nam,targString));
    bouton(i) = sum(regexp(nam,'bouton'));
    glom(i) = sum(regexp(nam,'glom'));
    center(i) = sum(regexp(nam,'center'));

    

    
end

ofID = frag | part | isCell | glom | center;

%% place in structure
nameProps.names = objNames;
nameProps.cellNum = cellNum;
nameProps.allIDs = allIDs;
nameProps.tcr = tcr>0;
nameProps.lin = lin>0;
nameProps.rgc = (cellNum>=1000 )& (cellNum<10000);
nameProps.postNum = postNum;
nameProps.cell = isCell;
nameProps.part = part;
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


%% Get edges
%edges
numEdge = 0;
numJuncs = 0;
edges = [];
juncs = [];
for i = 1:length(nameProps.names)
    if nameProps.part(i)
       parts = nameProps.allIDs{i};
       if length(parts) == 2
                numEdge = numEdge+1;
  
       edges(numEdge,:) = [parts i];
       end
    end
     if nameProps.junc(i)
       parts = nameProps.allIDs{i};
       if length(parts) == 2
               numJuncs = numJuncs+1;
   
       juncs(numJuncs,:) = [parts i];
       end
     end
end
% 
% postList = sort(unique(edges(:,1)));
% preList = sort(unique(edges(:,2)));

nameProps.edges = edges;
nameProps.juncs = juncs;

%% Make connectivity matrix.
if isempty(edges)
    preID = [];
    postID = [];
else
preID = unique(edges(:,1));
postID = unique(edges(:,2));
end





function[nep,syn] = swc2nep(fileName)


formatSpec = '%f%f%f%f%f%f%f%[^\n]';
delimiter = {',',' '};
%headerRows = 11;


%% count header lines
for i = 1:100;
    fid = fopen(fileName);
    sDat = textscan(fid,formatSpec, Inf, 'Delimiter', delimiter, 'HeaderLines', i);
    fclose(fid);
    if length(sDat{1})
        headerRows = i;
        break
    end
end


%% swc file to sDat
goodRead = 1;
try
    fid = fopen(fileName);
    sDat = textscan(fid,formatSpec, Inf, 'Delimiter', delimiter, 'HeaderLines', headerRows);
    fclose(fid);
catch err
    goodRead = 0;
    disp('problem reading swc file');
end


%% sDat to nodes and edges


nodeType = sDat{2};

isPre = find(nodeType == 7);
isPost = find(nodeType == 8);
isSyn = cat(1,isPre, isPost);
isSkel = setdiff([1:length(nodeType)],isSyn);
isLinked = setdiff(isSkel,1);
newNodes = zeros(length(nodeType),1);
newNodes(isSkel) = 1:length(isSkel);

if goodRead

    nep.edges = cat(2,sDat{1}(isLinked),sDat{7}(isLinked));
    useNodes = unique(nep.edges);
    newNodes = zeros(length(nodeType),1);
    newNodes(useNodes) = 1:length(useNodes);
    nep.edges = newNodes(nep.edges);
    
    nep.pos = cat(2,sDat{3}(useNodes),sDat{4}(useNodes),sDat{5}(useNodes));
    nep.nodeRad = sDat{6}(useNodes);
    nep.nodes = sDat{1}(useNodes);
    nep.edgeRad = (nep.nodeRad(nep.edges(:,1)) + nep.nodeRad(nep.edges(:,2)))/2;
else
    nep = [];
end


syn.pos = cat(2,sDat{3}(isSyn),sDat{4}(isSyn),sDat{5}(isSyn));
syn.pre = zeros(length(isSyn),1);
syn.pre(1:length(isPre)) = 1;
syn.post = zeros(length(isSyn),1);
syn.post(length(isPre)+1:length(isSyn)) = 1;
syn.nearest = sDat{7}(isSyn);




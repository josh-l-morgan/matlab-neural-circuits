
global glob tis
app = glob.app;



view(glob.ax,0,270)


cids = [1 10];


for c = 1:length(cids)
    cid = cids(c);
    fileName = sprintf('%s%d.mat',glob.useFvDir,cid);
    fv = loadFV(fileName);


    fileNameLab = sprintf('%s%d_labels_44.mat',glob.useFvDir,cid);
    data = load(fileNameLab);
    
    fv_1_vertices = data.fv_1_vertices;
    vertLab = fv_1_vertices(:,4);
    vertLab(vertLab==0) = 3;
    cMap = [0 0 0; 1 0 1; 0 0 0];

    newG = length(glob.g)+1;
    idx = find(tis.cids==cid);
    %% make new group
    glob.g(newG) = glob.defaultG;
    glob.g(newG).idx = idx;
    glob.g(newG).cid = cid;
    glob.g(newG).col = [1 1 1];
    glob.g(newG).alph = 1;
    glob.g(newG).show = 1;
    glob.g(newG).name = sprintf('group %d',cid);


    %% Show groups
    showCellGroup(app, newG)
    updateGroups(app)
    %edit2_function(app)

    pIdx =  glob.g(newG(1)).patch;
    set(pIdx,"FaceVertexCData",cMap(vertLab,:))
    set(pIdx,"FaceColor",'interp')
    set(pIdx,"EdgeColor",'none')

end

%%
use_on_group = [6, 7, 8, 9, 10, 11, 12, 14, 17, 20, 21, 23];
use_off_group = [1, 2, 3, 4, 5, 15, 18, 19];


% break tis.syn

cid = 1;

all_nodes = fv_1_vertices;
prox_nodes = all_nodes(all_nodes(:, 4) == 1, 1:3);
prox_reordered1 = prox_nodes(:, [2, 3, 1]);
prox_reordered1=prox_reordered1- [11,-10,0];

dist_nodes =  all_nodes(all_nodes(:, 4) == 2, 1:3);
dist_reordered1 = dist_nodes(:, [2, 3, 1]);
dist_reordered1=dist_reordered1- [11,-10,0];

all_nodes = fv_10_vertices;
prox_nodes = all_nodes(all_nodes(:, 4) == 1, 1:3);
prox_reordered10 = prox_nodes(:, [2, 3, 1]);
prox_reordered10=prox_reordered10- [11,-10,0];

dist_nodes =  all_nodes(all_nodes(:, 4) == 2, 1:3);
dist_reordered10 = dist_nodes(:, [2, 3, 1]);
dist_reordered10 = dist_reordered10 - [11,-10,0];


prox_syn_ind = [];
dist_syn_ind = [];

for i = 1:length(tis.syn.pre)
    if tis.syn.post(i) == 1
        if tis.syn.preClass(i) == 7
            r = i;
            pos_syn = tis.syn.pos(i,:);
            pos_syn = repmat(pos_syn, size(prox_reordered1, 1), 1);
            % Calculate the Euclidean distance between corresponding rows
            diff = sqrt(sum((prox_reordered1 - pos_syn).^2, 2));
            min_dif_prox = min(diff);
            pos_syn = tis.syn.pos(i,:);
            pos_syn = repmat(pos_syn, size(dist_reordered1, 1), 1);
            diff = sqrt(sum((dist_reordered1 - pos_syn).^2, 2));
            min_dif_dist = min(diff);
            if min_dif_prox < min_dif_dist
                prox_syn_ind = [prox_syn_ind i];
                min_dif_prox
            end
            if min_dif_prox > min_dif_dist
                dist_syn_ind = [dist_syn_ind i];
                min_dif_dist
            end   
        end
    end
    
    if tis.syn.pre(i) == 1
        r = i;
        pos_syn = tis.syn.pos(i,:);
        pos_syn = repmat(pos_syn, size(prox_reordered1, 1), 1);
        % Calculate the Euclidean distance between corresponding rows
        diff = sqrt(sum((prox_reordered1 - pos_syn).^2, 2));
        min_dif_prox = min(diff);
        pos_syn = tis.syn.pos(i,:);
        pos_syn = repmat(pos_syn, size(dist_reordered1, 1), 1);
        diff = sqrt(sum((dist_reordered1 - pos_syn).^2, 2));
        min_dif_dist = min(diff);
        if min_dif_prox < min_dif_dist
            prox_syn_ind = [prox_syn_ind i];
            min_dif_prox
        end
        if min_dif_prox > min_dif_dist
            dist_syn_ind = [dist_syn_ind i];
            min_dif_dist
        end   
    end
    
    if tis.syn.post(i) == 10
        if tis.syn.preClass(i) == 7
            r = i;
            pos_syn = tis.syn.pos(i,:);
            pos_syn = repmat(pos_syn, size(prox_reordered10, 1), 1);
            % Calculate the Euclidean distance between corresponding rows
            diff = sqrt(sum((prox_reordered10 - pos_syn).^2, 2));
            min_dif_prox = min(diff);
            pos_syn = tis.syn.pos(i,:);
            pos_syn = repmat(pos_syn, size(dist_reordered10, 1), 1);
            diff = sqrt(sum((dist_reordered10 - pos_syn).^2, 2));
            min_dif_dist = min(diff);
            if min_dif_prox < min_dif_dist
                prox_syn_ind = [prox_syn_ind i];
            end
            if min_dif_prox > min_dif_dist
                dist_syn_ind = [dist_syn_ind i];
            end   
        end
    end
    
    
if tis.syn.pre(i) == 10
        r = i;
        pos_syn = tis.syn.pos(i,:);
        pos_syn = repmat(pos_syn, size(prox_reordered10, 1), 1);
        % Calculate the Euclidean distance between corresponding rows
        diff = sqrt(sum((prox_reordered10 - pos_syn).^2, 2));
        min_dif_prox = min(diff);
        pos_syn = tis.syn.pos(i,:);
        pos_syn = repmat(pos_syn, size(dist_reordered10, 1), 1);
        diff = sqrt(sum((dist_reordered10 - pos_syn).^2, 2));
        min_dif_dist = min(diff);
        if min_dif_prox < min_dif_dist
            prox_syn_ind = [prox_syn_ind i];
        end
        if min_dif_prox > min_dif_dist
            dist_syn_ind = [dist_syn_ind i];
        end   
    end
    
end





edges = [];
pre = [];
post = [];
obID = [];
synType = [];
preClass = [];
postClass = [];
pos = [];
synPosRaw = [];
synPosDS = [];

in = 1;
% Loop through the elements
for j = 1:length(tis.syn.pre)
    if ~ismember(j, dist_syn_ind)
        edges(in, :) = tis.syn.edges(j, :);
        pre(in, :) = tis.syn.pre(j, :);
        post(in, :) = tis.syn.post(j, :);
        obID(in, :) = tis.syn.obID(j, :);
        synType(in, :) = tis.syn.synType(j, :);
        preClass(in, :) = tis.syn.preClass(j, :);
        postClass(in, :) = tis.syn.postClass(j, :);
        pos(in, :) = tis.syn.pos(j, :);
        synPosRaw(in, :) = tis.syn.synPosRaw(j, :);
        synPosDS(in, :) = tis.syn.synPosDS(j, :);
        in = in + 1;
    end
end

% save the above variables in syn

syn = struct();

syn.edges = edges;
syn.pre = pre;
syn.post = post;
syn.obID = obID;
syn.synType = synType;
syn.preClass = preClass;
syn.postClass = postClass;
syn.pos = pos;
syn.synPosRaw = synPosRaw;
syn.synPosDS = synPosDS;
syn.synProp = tis.syn.synProp ; 
syn.order = tis.syn.order ;
syn.prePostTypeName = tis.syn.prePostTypeName ;
syn.prePostTypeRule = tis.syn.prePostTypeRule;

obI = tis.obI;
cids = tis.cids;
cells = tis.cells;
dat = tis.dat;
report = tis.report;


tis = struct();

tis.obI = obI;
tis.cids = cids;
tis.cells = cells;
tis.syn = syn;
tis.dat = dat;
tis.report = report;


% Define the folder and file name
folder = 'Y:/Active/morganLab/pratyushsRetina/jxQ_cellnav_figure/Volumes/JxQ_vol2_highMed/Analysis/';
filename = 'tis_prox.mat';

% Construct the full file path
fullFilePath = fullfile(folder, filename);

% Save the variable 'tis' to the specified folder
save(fullFilePath, 'tis');








 
        
    
    




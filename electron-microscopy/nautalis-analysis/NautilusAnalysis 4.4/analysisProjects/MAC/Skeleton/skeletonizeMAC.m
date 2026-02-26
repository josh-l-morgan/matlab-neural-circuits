
MPN = 'Z:\joshm\LGNs1\Exports\MAC_export\MAC_merge_mat\';

load([MPN 'obI.mat'])
load([MPN 'dsObj.mat'])

%% Select cells
allEdges = obI.nameProps.edges;
allEdges(allEdges == 11007) = 1007;
preTarg = preTo(allEdges,1007);
postTarg = postTo(allEdges,1007);
%%Reference lists
showCell = [1007 11007];

% isPost = preTarg(:,1);
% isPre = postTarg(:,1);

allPost = unique(allEdges(:,1));
allPost = allPost(allPost>0);
allPre = unique(allEdges(:,2));

cells = obI.cell.name(obI.cell.isCell>0);
cellList = allPost;
% cellList = intersect(cells,allPost);
% cellList = 2011;

disp(sprintf('Cell List = %s',num2str(cellList)));
TPN = [MPN 'skel\'];
TPNview = [TPN 'view\'];
TPNmat = [TPN 'mat\'];
if ~exist(TPN,'dir'),mkdir(TPN);end
if ~exist(TPNview,'dir'),mkdir(TPNview);end
if ~exist(TPNmat,'dir'),mkdir(TPNmat);end

for i = 1:length(cellList)
    clf
    cellTarg = cellList(i);
    disp(sprintf('skeletonizing cell %d (%d of %d)',cellTarg,i,length(cellList)))
    
    skelFile = sprintf('%s%d.mat',TPNmat,cellTarg);
    viewFile = sprintf('%s%d.png',TPNview,cellTarg);
    if 1;%~exist(skelFile,'file')
        
        startTime = clock;
        
        rawObjectSubs = getCellSubs(obI,dsObj,cellTarg);
        rawSeed   = ceil(getSeed(obI,cellTarg));
        fixSeed = rawSeed.* obI.em.res .* [8 8 1]  ./obI.em.dsRes/1000;
        dsamp = [1 1 1];
        if (size(rawObjectSubs,1)>0) 
            pass = 1;
            cellFailed = [];
            
            objectSubs = double(downSampSub(rawObjectSubs,dsamp));
            %probableDownSamp = rawSeed./(median(rawObjectSubs)./[4 4 4])
             
            if isempty(rawSeed)
                seedSub = median(objectSubs,1);
            else
                seedSub = round(fixSeed./[4 4 4]);
            end
        
            
            cellStruct = subs2arbor(objectSubs,seedSub);
            cellStruct.voxSize = obI.em.dsRes.*dsamp;
            
%             try
%                 
%             catch err
%                 err
%                 pass = 0;
%                 cellFailed = [cellFailed cellTarg]
%             end
            
            if pass
                cellView = cellStruct.sideViews{1};
                skel = cellStruct.skel;
                
                
                save(skelFile,'cellStruct')
                imwrite(cellView,viewFile)
                image(cellView)
%                 
%                 cellArbors.cellName(i) = cellTarg;
%                 cellArbors.arbors(i) = rmfield(cellStruct.arbor,'vox');
            end
            
            
            stopTime = clock;
        
        clf
   % plotSkeletons(cellArbors.arbors(i),seedSub);
%     fileName = sprintf('%smorphOri_%03.0f.png',morphDir,cellTarg)
%     saveas(gcf,fileName,'png')
        
        end
        
    end
    
    
end



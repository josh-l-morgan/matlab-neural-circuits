


%% Skeletonize with shortest paths
cellList = obI.cell.name;
cellList = intersect(obI.cell.name,[1:999]);

disp(sprintf('Cell List = %s',num2str(cellList)));
TPN = [MPN 'skel\'];
TPNview = [TPN 'view\'];
TPNmat = [TPN 'mat\'];
if ~exist(TPN,'dir'),mkdir(TPN);end
if ~exist(TPNview,'dir'),mkdir(TPNview);end
if ~exist(TPNmat,'dir'),mkdir(TPNmat);end
morphDir = [MPN 'morph\'];
if ~exist(morphDir,'dir'),mkdir(morphDir),end


for i = 1:length(cellList)
    clf
    cellTarg = cellList(i);
    disp(sprintf('skeletonizing cell %d (%d of %d)',cellTarg,i,length(cellList)))
    
    skelFile = sprintf('%s%d.mat',TPNmat,cellTarg);
    viewFile = sprintf('%s%d.png',TPNview,cellTarg);
    if ~exist('skelFile','file')
        
        startTime = clock;
        
        rawObjectSubs = getCellSubs(obI,dsObj,cellTarg);
        rawSeed   = ceil(getSeed(obI,cellTarg));
        
        if (size(rawObjectSubs,1)>1000) & (~isempty(rawSeed))
            pass = 1;
            cellFailed = [];
            
            objectSubs = double(downSampSub(rawObjectSubs,[4 4 4]));
            probableDownSamp = rawSeed./(median(rawObjectSubs)./[4 4 4])
            seedSub = round(rawSeed./[32 32 16]);
            cellStruct = subs2arbor(objectSubs,seedSub);
            
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
                
                cellArbors.cellName(i) = cellTarg;
                cellArbors.arbors(i) = rmfield(cellStruct.arbor,'vox');
            end
            
            
            stopTime = clock;
        
        clf
    plotSkeletons(cellArbors.arbors(i));
    fileName = sprintf('%smorphOri_%03.0f.png',morphDir,cellTarg)
    saveas(gcf,fileName,'png')
        
        end
         
        
        
    end
    
   
    
    
end

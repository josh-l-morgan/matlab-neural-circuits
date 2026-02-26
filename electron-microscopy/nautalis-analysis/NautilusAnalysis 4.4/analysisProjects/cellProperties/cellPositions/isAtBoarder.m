function[useNames] = isCellBody();




load('MPN.mat')
if ~exist('MPN','var')
    MPN = GetMyDir;
end

synDir = [MPN 'synPos3\'];
if ~(exist(synDir,'dir')),mkdir(synDir); end

if ~exist('obI','var') | ~exist('dsObj','var')
    disp('loading')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
end


%
% mot = getMotifs(obI);
% syn = getSynMat(obI);


%%
clear vol
anchors  = obI.cell.anchors;

dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;
anchors(:,1) = anchors(:,1)*dSamp(1);
anchors(:,2) = anchors(:,2)*dSamp(2);
anchors(:,3) = anchors(:,3)*dSamp(3);

r = 5;
maxVol = 4/3 * pi * r^4;
clear mins maxes


for o = 1:length(obI.cell.name)
    
    ids = obI.cell.obIDs{o};
    subs = [];
    for i = 1:length(ids)
        sub = dsObj(ids(i)).subs;
        subs = cat(1,subs,sub);
    end
    
    maxS = max(subs,[],1);
     if isempty(subs)
             mins(i,:) = [0 0 0];
            maxes(i,:) = [inf inf inf];
        else        
         
        inds = sub2ind(maxS,subs(:,1),subs(:,2),subs(:,3));
        inds = unique(inds);
        [y x z] = ind2sub(maxS,inds);
        subs = [y x z];
        
        anc = anchors(o,:);
        
        scatter(subs(:,1),subs(:,2),'.')
        
        hold on 
        
        scatter(anc(1),anc(2),'o')
        hold off
       

            
           mins(o,:) = min(subs,[],1);
           maxes(o,:) = max(subs,[],1);
        end
   
    
    %     b = hist(dist,[1:100])
    %     bar(b(1:end-1));
    %     ylim([0 50000])
    %     pause
    %     sum(dist<5)
    %
    
end


useNames = obI.cell.name(vol>5);



























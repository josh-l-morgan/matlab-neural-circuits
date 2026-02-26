
clear all

MPN = 'Z:\joshm\LGNs1\Exports\MAC_export\MAC_merge_mat\';
DPN = [MPN 'data\'];

load([DPN 'topoPairDists4.mat'])
load([MPN 'obI.mat']);



goodComp = 0;
clear eucDists pLengths compileLengths
minSpac = .1;
for t = 1:length(topoPairDists)
    pLength = topoPairDists(t).pLength;
    eucDist = topoPairDists(t).eucDist;
    sPos = topoPairDists(t).sPos;
    
    hold off
    scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'o','filled','k')
    hold on
    
    for p = 1:length(pLength)
        path = topoPairDists(t).keepPaths{p};
        %path = flipud(path);
        lengths = topoPairDists(t).keepLengths{p};
        if ~isempty(lengths)
            
           hold off
            scatter3(sPos(:,1),sPos(:,2),sPos(:,3),'.','k');
            hold on
            
            L = 0;
            c = 0;
            lastN = path(1,1);
            clear shortPath shortLengths
            
            
            for n = 1:size(path,1)
%                 plot(sPos(path(n,:),1)', sPos(path(n,:),2)','r');
%                 pause(.01)
                dist = sqrt((sPos(lastN,1)-sPos(path(n,2),1)).^2 + ...
                    (sPos(lastN,2)-sPos(path(n,2),2)).^2 + ...
                    (sPos(lastN,3)-sPos(path(n,2),3)).^2);
                 
%                 plot3(sPos(path(n,:),1),sPos(path(n,:),2),...
%                     sPos(path(n,:),3),'g','linewidth',4)
                
                if (dist > minSpac) | (n == size(path,1))
                    c = c+1;
                    shortLengths(c) = dist;
                    shortPath(c,:) = [lastN path(n,2)];
                    lastN = path(n,2);
                    
                end
            end
                            pause(.01)

            %shortLengths,pause
            
    plot3(sPos(shortPath(:,1),1),sPos(shortPath(:,1),2),...
                   sPos(shortPath(:,1),3),'r','linewidth',2)
   plot3(sPos(shortPath(end,:),1),sPos(shortPath(end,:),2),...
                   sPos(shortPath(end,:),3),'g','linewidth',2)
               
          scatter3(sPos(shortPath(1,1),1),sPos(shortPath(1,1),2),...
                   sPos(shortPath(1,1),3),'g','o','filled')
     scatter3(sPos(shortPath(1,1),1),sPos(shortPath(1,1),2),...
                   sPos(shortPath(1,1),3),'b','o','filled')            
            pause(.01)

           
             'lengths'
            pause(.1)
            goodComp = goodComp+1
            compileLengths(goodComp) = sum(shortLengths);
            eucDists(goodComp) = eucDist(p);
            fromCell(goodComp) = topoPairDists(t).cellName;
            
            
        end
        
    end
    
end

hold off
scatter(eucDists,compileLengths,'k','o','filled')
hold on
plot([0 30],[0 30],'r')
hold off



usedCells = unique(fromCell);
histCells = hist(fromCell,usedCells);










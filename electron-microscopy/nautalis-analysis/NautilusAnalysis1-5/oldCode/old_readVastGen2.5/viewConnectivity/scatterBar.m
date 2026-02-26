
function[cellProp scaleProp] = scatterBar(rawProp,anchors,barHeight);
%%
%cellProp = (PostConToSeed(:,1:2))
cellProp = rawProp;
barPos = anchors;

[cellNum propNum] = size(cellProp);

barWidth = 4;
barSpace = 1;
if ~exist('barHeight','var')
    barHeight = 50;
end

barCol = [.5 0 0 ; 0 .5 0; 0 0 .3];
backThick = .3;

scaleProp = barHeight/max(cellProp(:));
cellProp = cellProp*scaleProp;
cellProp(:,2) = cellProp(:,2) ;
shiftY = (cellProp(:,1)-cellProp(:,2))/2 ;
barPos(:,1) = barPos(:,1) - shiftY;

showBack = 1;


clf
hold on
for p = 1:cellNum
    
               


    for b = 1:propNum;
            if b>1;
                flip = -1;
            else
                flip = 1;
            end
                plotBackY = [barPos(p,1) barPos(p,1) + cellProp(p,b) * flip + backThick * flip];
                plotBackX = barPos(p,2);
               
                if showBack
                plot([plotBackX plotBackX], plotBackY,'LineWidth',barWidth+backThick*4,'color','w');
                end
                
            plotY = [barPos(p,1) barPos(p,1) + cellProp(p,b)* flip];
            plotX = barPos(p,2) ;
            if  cellProp(p,b)
                plot([plotX plotX], plotY,'LineWidth',barWidth,'color',barCol(b,:));
            end
            
            %plot(plotX, plotY,'LineWidth',1,'color',axCol);
    end
    
end


% 
% scatter(anchors(:,2),anchors(:,1),2,'k','filled','o',...
%     'MarkerEdgeColor','w','LineWidth',1)

            pause(.01)
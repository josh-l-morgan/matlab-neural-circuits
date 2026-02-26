function[shuffleGraph,newOrder] = reshuffleSimilarity(simGraph);


rateGraph = simGraph * 0 + 1;
[y x] = find(rateGraph);
mid = (y+x)/2;
dists = sqrt((y -mid).^2 + (x - mid).^2);
rateGraph(:) = dists;
image(rateGraph*3);

shuffleGraph = simGraph;
postNum = size(shuffleGraph,1);

newOrder = 1:postNum;

image(shuffleGraph * 10)

reps = 10000000;
maxNoChange = 10000;
countNoChange = 0;
preRand = 0;
for rep = 1:reps
    
    swap1 = fix(rand*postNum)+1;
    swap2 = fix(rand*postNum)+1;
    
    newGraph = shuffleGraph;
    newGraph(:,swap1) = shuffleGraph(:,swap2);
    newGraph(:,swap2) = shuffleGraph(:,swap1);
    
    tempLine = newGraph(swap1,:);
    newGraph(swap1,:) = newGraph(swap2,:);
    newGraph(swap2,:) = tempLine;
    
    newRating = (newGraph .* rateGraph);
    oldRating = (shuffleGraph .* rateGraph);
    if(rep<preRand)
        countNoChange = 0;
        shuffleGraph = newGraph;
        tempOrder = newOrder(swap1);
        newOrder(swap1) = newOrder(swap2);
        newOrder(swap2) = tempOrder;
        
    else
    if (sum(newRating(:)) < sum(oldRating(:)))
        countNoChange = 0;
        shuffleGraph = newGraph;
        tempOrder = newOrder(swap1);
        newOrder(swap1) = newOrder(swap2);
        newOrder(swap2) = tempOrder;
        
    else
        countNoChange = countNoChange+1;
    end
    end
    
    if ~mod(rep,100)
        image(shuffleGraph*250/max(shuffleGraph(:))),pause(.01)
        end
    
    
    if countNoChange > maxNoChange
        'no change'
        break
    end
end

    
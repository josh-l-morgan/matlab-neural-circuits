
cent = [12805, 19582, 3860]

pos = ...
'(14924, 20813, 4038)'
res = [.016 .016 .03];


startP = regexp(pos,'(');
stopP = regexp(pos,')');
clear pointName allPoints pointDist
for i = 1: length(startP)
    pointName{i} = pos(startP(i)+1:stopP(i)-1)
    pointPos= str2num(pointName{i})
   allPoints{i}= pointPos
   pointDist(i) = sqrt(((pointPos(1)-cent(1))*res(1)).^2 + ((pointPos(2)-cent(2))*res(2)).^2 + ((pointPos(3)-cent(3))*res(3)).^2 )
end


pointDist




%%


cent = [13058, 19560, 3872]

%      pos = {}

res = [.016 .016 .03];
clear allPoints pointDist
pos = {};
allPoints = {};
pointDist = {}
for i = 1: size(pos,1)
    pointName = pos{i};
    if isempty(pointName)
         pointDist(i,1) = NaN
        
        
    else
        pointPos = str2num(pointName(2:end-1));
       allPoints{i}= pointPos
   pointDist(i,1) = sqrt(((pointPos(1)-cent(1))*res(1)).^2 + ((pointPos(2)-cent(2))*res(2)).^2 + ((pointPos(3)-cent(3))*res(3)).^2 )

    
    end
  end



%%


cent = [13058, 19560, 3872]

%      pos = {}

res = [.016 .016 .03];
clear allPoints pointDist
allPoints = {};
pointDist = []
for i = 1: size(pos,2)
    pointName = pos{i};
    if isempty(pointName)
         pointDist(1,i) = NaN
        
        
    else
        pointPos = str2num(pointName(2:end-1));
       allPoints{i}= pointPos
   pointDist(1,i) = sqrt(((pointPos(1)-cent(1))*res(1)).^2 + ((pointPos(2)-cent(2))*res(2)).^2 + ((pointPos(3)-cent(3))*res(3)).^2 )

    
    end
  end

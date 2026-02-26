
cent = [11727, 19706, 4469]

pos = ...
'(11727, 19706, 4469)	(11278, 19795, 4834)	(11158, 19761, 4944)	(11718, 19940, 5167)	(11731, 20187, 5006)	(11833, 20521, 4414)	(11745, 20467, 4301)'
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










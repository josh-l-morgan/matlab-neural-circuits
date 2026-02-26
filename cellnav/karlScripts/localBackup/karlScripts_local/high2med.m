function outCoords=high2med(X,Y,Z,inv)
inv
if ~exist('inv','var')
    vert=0;
else
    if strcmp(inv,'inv')
        vert=1;
    else
        vert=0;
    end
end
tfm=affine2d();
tfm.T=[4.70655613870859,0.192569634745259,0;-0.244283296195143,4.89254876488754,0;-62975.9722844041,-59712.5119665449,1]
if vert==0
[x,y]=tfm.transformPointsInverse(X,Y);
else
[x,y]=tfm.transformPointsForward(X,Y);    
end
rawCoords=[x,y,Z];
outCoords=uint16(rawCoords)

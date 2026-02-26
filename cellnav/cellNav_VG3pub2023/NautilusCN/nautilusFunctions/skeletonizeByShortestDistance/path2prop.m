function[squeezedPath] = path2prop(surfVox,prop,reps)

prop = prop(:)';

%% for each vox on path from tip, enter distance of longest source tip

path2pathOld = surfVox.seedPath;
%path2path = passMaxAndDist2(surfVox,prop,reps);
path2path = conMat2allShortestProp(surfVox,prop,reps);
subs = surfVox.subs;

%combine pred lists
newPred = path2path.pred;
notChanged = find(newPred<1);
newPred(notChanged) =  path2pathOld.pred(notChanged);
predLength = path2path.predLength;
predLength(notChanged) = path2pathOld.predLength(notChanged);

%showMaxProp(moveObj,notChanged.*seedPath.dist'/5+20);

%Get new path
[pathLN] = drawPathDist(newPred,path2pathOld.dist,predLength);

[sP sIDX] = sort(pathLN.lengths,'descend');
%showPathL3D3(surfVox.subs,pathLN,pathLN.voxLengths,sIDX(1:700))



squeezedPath = pathLN;
squeezedPath.pred = newPred;
squeezedPath.seed = surfVox.firstSeed;
squeezedPath.prop = path2path.prop;
squeezedPath.path2path = path2path;















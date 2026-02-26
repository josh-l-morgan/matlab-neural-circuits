function[col] = showPath3D(vox,path)


%%
if ~exist('projectDim','var')
    projectDim = 1;
end

dims = [1:3];
dims = dims(dims ~= projectDim);

subs = vox.subs;
vsubs = vox.subs;
offset = min(subs,[],1)-1;

if ~exist('dims','var')
    dims = [1 2];
end

for i = 1:3
    subs(:,i) = subs(:,i) - offset(i) + 1;
    vsubs(:,i) = vsubs(:,i) - offset(i) + 1;
end

%% Render 3D
hold off
clf
downSamp = 4;
renderProps.smooth = 0;
renderProps.resize = 1;
renderProps.smoothPatch = 0;

smallSub = shrinkSub(vsubs,downSamp);
%smallSub = smallSub(:,flipDim);

fv = subVolFV(smallSub,[],renderProps);
%fv.vertices = fv.vertices * downSamp;
fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;

[p] = renderFV(fv,[0 0 1],.4);
view([0 0])
axis off
pause(.01)
hold on



for b = 1:length(path.bones)
    
    runE = subs(path.bones(b).nodes,:);
    scatter3(runE(:,1),runE(:,2),runE(:,3),300,'r','.')
   
    plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','g')
    
end

hold off



function[fv] = ob2OBJ(obI,dsObj,showCellNames)



%%
uCellId = 201; i = 1;
targ = find(obI.cell.name==uCellId(i));
obTarg = obI.cell.obIDs{targ};
obSub = [];

for o = 1:length(obTarg)
    if obTarg(o)<=length(dsObj)
        sub = double(dsObj(obTarg(o)).subs);
        if ~isempty(sub)
            obSub = cat(1,obSub, sub);
        end
    end
end
    
%make vol
obSub = round(obSub);
obSub(:,1) = obSub(:,1)- min(obSub(:,1))+1;
obSub(:,2) = obSub(:,2)- min(obSub(:,2))+1;
obSub(:,3) = obSub(:,3)- min(obSub(:,3))+1;


obSiz = max(obSub,[],1)+1;

obInd = sub2ind(obSiz,obSub(:,1),obSub(:,2),obSub(:,3));

obVol = zeros(obSiz);
obVol(obInd) = 1;


%%
        

%%
              %[f,v]=voxel_bnd_faces(obSub,[1 1 1],[0 0 0],1);

    
              %%
tic
fv=isosurface(obVol,0.5);
toc
% 
% tic
% fv = obConFace(obSub)
% toc


[f v] = isosurface(obVol,0.5);

clf

p = patch(fv);
view(30,-15);
axis vis3d;
colormap copper
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud
        

%%


%%
    filename='D:\LGNs1\Segmentation\VAST\S4\joshm\exports\test\testSmall.obj';
    
              objectname='hi';
              disp(['Saving ' filename ' as Wavefront OBJ.....']);
              vertface2obj(v,f,filename,objectname);
vertface2obj




        
        
        
function[fv] = subVolFVproject(obSub,saveName,renderOb)



%% make vol
obSub = round(obSub);
% obSub(:,1) = obSub(:,1)- min(obSub(:,1))+1;
% obSub(:,2) = obSub(:,2)- min(obSub(:,2))+1;
% obSub(:,3) = obSub(:,3)- min(obSub(:,3))+1;


obSiz = max(obSub,[],1)+1;

obInd = sub2ind(obSiz,obSub(:,1),obSub(:,2),obSub(:,3));

obVol = zeros(obSiz);
obVol(obInd) = 1;

%% project
for i = 1:size(obSub,1)
   obVol(obSub(i,1),obSub(i,2),1:obSub(i,3)) = 1;     
end

%%
     
fv=isosurface(obVol,0.5);
if ~exist('renderOb','var')
    renderOb = 0;
end

 if renderOb 

p = patch(fv);
view(30,-15);
axis vis3d;
colormap copper
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud
 end  

%%


%%
if exist('fileName','var')
              objectname='hi';
              disp(['Saving ' filename ' as Wavefront OBJ.....']);
              vertface2obj(v,f,filename,objectname);
end


        
        
        
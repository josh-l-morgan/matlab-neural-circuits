


%clear all

KPN=GetMyDir
Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));
yxum=0.103;
zum=0.3;
%load('cmap.mat')
%set(gcf,'Colormap',cmap)


for k = 1:1%size(Kdir,1)
    k
    clear Branch 
    TPN = [KPN '\' Kdir(k).name '\'];
    if exist([TPN 'BranchS.mat'])
        load([TPN 'BranchS.mat'])
        load([TPN 'Use.mat'])
        Mids=Use.Mids;
        BranchMap=Branch.BranchMap;
            
           DrawSeg=fix(Mids)+1;
           DrawBranch=zeros(max(DrawSeg(:,1)),max(DrawSeg(:,2)));
           UseMid=find(BranchMap);
           for ds=1:size(UseMid,1)
                DrawBranch(DrawSeg(UseMid(ds),1),DrawSeg(UseMid(ds),2))=50*Branch.DotCount(BranchMap(UseMid(ds)))+2;           
           end
            image(DrawBranch),pause(.1)
        
    
%         
%        if exist([TPN 'Cell.mat'])
%             load([TPN 'Cell.mat'])
%             cName=['P' Cell.Age '_' Cell.Type '_' Kdir(k).name '.tif']
%             Name=[KPN 'Images\DrawBranches\' cName];
%         else
%             Name=[KPN 'Images\DrawBranches\' Kdir(k).name '.tif'];
%         end
%         if isempty(find(Name=='?'));
%             imwrite(DrawBranch,gaymap,Name,'Compression','none')
%         end
      
    
    
    end
end
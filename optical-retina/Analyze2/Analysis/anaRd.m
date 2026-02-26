function[]=anaDF(TPN,DPN)


load([TPN 'Dots.mat'])



%% Find appropriate mean faces

TSphere=ones(11,11,11);
[ty tx tz]=find3(TSphere);
TSv=[ty tx tz];
Tdists=dist(TSv,[6 6 6]);
for d=1:100
   Near=TSv(Tdists<(d/10),:);
   Tvol(d)=size(Near,1);
   clear Conn Faces
   for n = 1:Tvol(d)
       Conn=dist(Near,Near(n,:));
       Faces(n)=6-sum(Conn==1);
   end
   meanFaces(d)=mean(Faces);
   histFaces(d,:)=hist(Faces,0:1:6);
end
c=0;
for v=1:max(Tvol)
   if ~isempty(find(Tvol==v))
       c=c+1;
       tvol(c)=v;
       v2f(c)=meanFaces(find(Tvol==v,1));
   end    
end
RoundFaces=interp1(tvol,v2f,1:max(tvol));



%% Run Dots
for i = 1 :Dots.Num
       
   Cent=Dots.Pos(i,:);
   Vox=Dots.Vox(i).Pos;   
   Dist = dist(Vox,Cent);
   MeanD= max(1,mean(Dist));
   DistN=Dist/MeanD;
   CentN=Cent/MeanD;
   VoxN=Vox/MeanD;
   
 %  MedN = median(DistN,1);
%   ExtN = max(DistN)/MedN;
   
   if size(Vox,1)>1, %if more then one voxel
       [Co, Sc,latent]=princomp(VoxN);
       Dots.Round.Var(i,:)=latent;
       Dots.Round.Long(i)=latent(1)/max(.1,latent(2));
       Dots.Round.Oblong(i)=max(abs(Sc(:,1)))/max(1,abs(max(Sc(:,2))));
       Dots.Round.SumVar(i)=sum(latent);
   else
       Dots.Round.Var(i,:)=[0;0;0];
       Dots.Round.Long(i)=0;
       Dots.Round.SumVar(i)=0;
       Dots.Round.Oblong(i)=1;
   end
   
  
   %% find surface area
   for v = 1:size(Dots.Vox(i).Pos,1)
       Conn=dist(Dots.Vox(i).Pos,Dots.Vox(i).Pos(v,:));
       Dots.Vox(i).Faces(v)=6-sum(Conn==1);
   end
   Dots.Round.histFaces(i,:)=hist(Dots.Vox(i).Faces,0:1:6);
   Dots.Round.meanFaces(i)=mean(Dots.Vox(i).Faces);
   Dots.Round.Compact(i)=Dots.Round.meanFaces(i)/RoundFaces(Dots.Vol(i));
   

end

save([TPN 'Dots.mat'],'Dots')


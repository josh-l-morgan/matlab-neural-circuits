

ID=uint8(zeros(Dots.ImSize));
P=find(pass')'; %% list of passing puncta
for i =1:length(P)
  ID(Dots.Vox(P(i)).Ind)=fix(rand*254)+1; 
end
imwriteNp([TPN 'pic\Ids'],ID,'ID')
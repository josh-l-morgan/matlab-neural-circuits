function[dS] = isoEM(S);

%% Down Sample
dsTest = imresize(S(:,:,1),1/8,'nearest');

[dys dxs] = size(dsTest);
dS = zeros(dys,dxs,size(S,3));
parfor i = 1:size(S,3)
   dS(:,:,i) = imresize(S(:,:,i),1/8,'nearest'); 
end

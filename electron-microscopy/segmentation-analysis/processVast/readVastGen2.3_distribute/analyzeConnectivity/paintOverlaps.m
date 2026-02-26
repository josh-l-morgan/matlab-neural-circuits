


clear all
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\'
SPN = [MPN 'manyTerSubs\'];
TPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\skel\'

load([SPN 'terSubs_Ds8_Ds1_Look0.mat'])


%%
% 
% downSamp = 1;
% lookDist = 1;
% 
% 
% 

%{
consider 1,1   1,3  2,3   10,1

3,3 - lots of apposition
4,1 - lots of appositions
8,2 - lots of appositions
10,8 - looks nice
12,1
13,1 lots of appo
13,7 lots of big appo
23,11
15,1 faciculation
%}
postList = ([108  129	109	117	162	131	116	137	130	135	106]);
preList = [1006        1009        1012        1014        1021        1025        1027        1028        1029        1030        1031        1032        1033 ...
            1034        1036        1037        1041        1050        1051        1053        1054        1055        1056];

        preList = preList(10);
        postList = postList(8);

allSubs = cat(1, terSubs.cells.terSubs);
dsMaxSub = max(allSubs,[],1)+10;
clear allSubs



viewVol = zeros(dsMaxSub,'uint8');
synMat = terSubs.synMat;
for pre = 1 : length(preList)
    
    disp(sprintf('running pre %d of %d',pre,length(preList)))
    diSub1 = terSubs.cells(find(terSubs.cellList == preList(pre))).terSubs;
    diInd1 = sub2ind(dsMaxSub,diSub1(:,1),diSub1(:,2),diSub1(:,3));
    
 for post = 1: length(postList)
        
       diSub2 = terSubs.cells(find(terSubs.cellList == postList(post))).terSubs;

        diInd2 = sub2ind(dsMaxSub,diSub2(:,1),diSub2(:,2),diSub2(:,3));
        
        viewVol = viewVol * 0;
        viewVol(diInd2) = 1;
        postSum = squeeze(sum(viewVol,1));
        viewVol = viewVol * 0;
        viewVol(diInd1) = 1;
        preSum = squeeze(sum(viewVol,1));

        
        
        
        sharedVox = intersect(diInd1,diInd2);
        viewVol = viewVol * 0;
        viewVol(sharedVox) = 1;
        overSum = squeeze(sum(viewVol,1));
        
        colSum = preSum* 250/max(preSum(:))+ (preSum>0) * 50;
        colSum(:,:,2) = overSum*50 + (overSum>0) * 50;
        colSum(:,:,3) = postSum* 250/max(postSum(:)) + (postSum>0)*50;
        
        colSum = preSum* 250/max(preSum(:))+ (preSum>0) * 50;
        colSum(:,:,2) = overSum*10 + (overSum>0) * 50;
        colSum(:,:,3) = postSum* 250/max(postSum(:)) + (postSum>0)*50;
        
        image(uint8(colSum));
        
        fileName = sprintf('appo_z%d_ax%d_post%d.png',dsMaxSub(3),preList(pre),postList(post));
        imwrite(uint8(colSum),[TPN fileName])
        synapseNumber = synMat(find(terSubs.preList == preList(pre)), find(terSubs.postList == postList(post)));
        disp(sprintf('pre %d, post %d, syn %d',preList(pre),postList(post),synapseNumber))

        
      
    end
end

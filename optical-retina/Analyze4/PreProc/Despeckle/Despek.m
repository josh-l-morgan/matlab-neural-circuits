 


Imax=max(I,[],4);
% subplot(2,1,1)
image((Imax))

[ys xs Chan zs]=size(I);

%% Despk
ch=1
Ch=zeros(ys+2,xs+2,zs+2);
Ch(2:ys+1,2:xs+1,2:zs+1)=I(:,:,ch,:);
Ones=ones(27,1);
sur=(1:27)~=14;
Ch1D=Ch*0;
for y = 2: ys+1, for x = 2 : xs+1, for z = 2 : zs+1
    Samp=Ch(y-1:y+1,x-1:x+1,z-1:z+1);
    Vox=Samp(14);
    Sur=Samp(sur);
    Val=min(max(Sur(:)),Vox);
    %Val=mean(Samp(:));
    ChD(y,x,z)=Val;   
end,end,y,end

image(max(Ch,[],3))
image(max(Ch1D,[],3))

%% Median filter
% 'median filtering',pause(.1)
% if kernelSize(1)>1
%     for i = 1: size(channel1,3)
%        channel1(:,:,i)=medfilt2(channel1(:,:,i),[kernelSize(1) kernelSize(1)]);
%     end
% end
% if kernelSize(2)>1
%     for i = 1: size(channel2,3)
%        channel2(:,:,i)=medfilt2(channel2(:,:,i),[kernelSize(2) kernelSize(2)]);
%     end
% end

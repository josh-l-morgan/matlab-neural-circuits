
%% image
subplot(3,1,1)
image(I(:,:,:,13))
image(I(:,:,2,14))


If=Shave(I,2);

subplot(3,1,2)
image(If(:,:,:,13))
image(If(:,:,2,14))

[Ifb, Ifm] =BandFind3D(If);

subplot(3,1,3);
image(Ifb(:,:,:,13))
image(Ifb(:,:,2,14)*3)


%%
Combo(:,:,1,:)=I(:,:,2,:);
Combo(:,:,2,:)=If(:,:,2,:);
Combo(:,:,3,:)=If(:,:,2,:)*0;

image(Combo(:,:,:,14)*3)
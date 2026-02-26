function drawimgplusgrid(imgnr,name,sgridx,sgridy,gridsquareangle,gridsquaresize)

image=imread(name); %load slice image
image=double(image)/255;
figure(imgnr);
imshow(image);
hold on;
plot(sgridx,sgridy,'k*');
plot(sgridx,sgridy,'r.');
plot(sgridx,sgridy,'wo');
plot(sgridx(2,3),sgridy(2,3),'g.');
plot(sgridx(2,2),sgridy(2,2),'b.');
plot(sgridx(1,2),sgridy(1,2),'m.');
plot(sgridx(3,2),sgridy(3,2),'y.');
plot(sgridx(4,4),sgridy(4,4),'k.');
plot(sgridx(5,5),sgridy(5,5),'c.');
draworientedsquares(sgridx,sgridy,gridsquareangle,gridsquaresize,'k');
hold off;
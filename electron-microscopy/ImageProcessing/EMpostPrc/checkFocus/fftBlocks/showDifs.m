for i = 1:size(difs,1);
    i


dif = difs(i,:);
subplot(1,2,1)
plot(dif/max(dif)),pause(.01)
%plot((dif-dif(1))/max(dif-dif(1))),pause(.01)
comp2to(i) = (dif(2)-dif(1))/(dif(end)-dif(1));
hold on
end
hold off
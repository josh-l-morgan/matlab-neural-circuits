function[sumObj] = showSurf(surfVox);

minMax = surfVox.minMax;

for i = 1:3
    moveObj(:,i) = surfVox.subs(:,i)-minMax{1}(i);
end

objVol = zeros(minMax{2},'uint8');
objVol(sub2ind(size(objVol),moveObj(:,1),moveObj(:,2),moveObj(:,3))) = 1;

for s = 1:3
    sumObj{s} = squeeze(sum(objVol,s));
    image(sumObj{s}*20), pause(.01)
end
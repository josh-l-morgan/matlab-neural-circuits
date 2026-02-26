function[pts2] = rigidPts(movePts,A)



alpha = ones(size(movePts,1),1);
tp1 = movePts(:,1) * A(1,1) + movePts(:,2) * A(1,2) + A(1,3) * alpha;
tp2 = movePts(:,1) * A(2,1) + movePts(:,2) * A(2,2) + A(2,3) * alpha;
pts2 = [tp1 tp2];
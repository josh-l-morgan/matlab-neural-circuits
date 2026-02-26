function[pts2] = rigidPtsStack(pts,As)


alpha = ones(size(pts,1),1);
tp1 = pts(:,1) .* As(pts(:,3),1,1) + pts(:,2) .* As(pts(:,3),1,2) + As(pts(:,3),1,3) .* alpha;
tp2 = pts(:,1) .* As(pts(:,3),2,1) + pts(:,2) .* As(pts(:,3),2,2) + As(pts(:,3),2,3) .* alpha;
pts2 = [tp1 tp2 pts(:,3)];







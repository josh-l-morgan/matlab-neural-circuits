function outpos=sortfit2d(inpos,fitpos)
%Assumes that inpos and fitpos are matrices of equal size (n,2).
%Returns a matrix which contains the rows of inpos re-ordered according to the minimal distances to points in fitpos.
%By Daniel Berger for MIT-BCS Seung, February 2010

nrpoints=size(inpos,1);

outpos=zeros(nrpoints,2);
used=zeros(nrpoints,1);
for i=1:1:nrpoints
  %compute distances to all points in fitpos
  dist=zeros(nrpoints,2);
  for j=1:1:nrpoints
    dx=fitpos(i,1)-inpos(j,1);
    dy=fitpos(i,2)-inpos(j,2);
    dist(j,1)=sqrt(dx*dx+dy*dy);
    dist(j,2)=j;
  end;
  dist=sortrows(dist,1);
  used(dist(1,2))=used(dist(1,2))+1;
  outpos(i,:)=inpos(dist(1,2),:);
end;

if (max(used)>1)
  disp('WARNING: sortfit2d did not result in a 1:1 fit!');
end;
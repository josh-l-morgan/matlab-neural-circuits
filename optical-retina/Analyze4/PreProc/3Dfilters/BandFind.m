function[IfMax,Maxs] = BandFind(Ip);
%%runs a range of gaussian filters and then maxes the results
%%Second output identifis the frequency of that max


[ys xs]=size(Ip);
Radius=20;

If=zeros(ys,xs,Radius,class(Ip));

for s = 1:Radius;
    siz=(s*2)-1;
   If(:,:,s)=imfilter(Ip,fspecial('gaussian',siz,siz/5)); 
end

% Maxi=max(If(:))
% for i = 1: Radius
%     i,image(If(:,:,i)*(100/Maxi)),pause(.01)
%     
% end

IfMax=max(If,[],3);
Same=uint16(Ip); Sames=Same; Maxs=Same;
for i = 1:Radius
   Same=uint16(If(:,:,i)==IfMax);
   Sames=Sames+Same;
   Maxs=Maxs+Same*i;
end
Maxs=Maxs./Sames;


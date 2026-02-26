function[subs]= wall1(subs,siz);


subs(subs<1) = 1;
for i = 1:length(siz)
    subs(subs(:,i)>siz(i),i) = siz(i);
end
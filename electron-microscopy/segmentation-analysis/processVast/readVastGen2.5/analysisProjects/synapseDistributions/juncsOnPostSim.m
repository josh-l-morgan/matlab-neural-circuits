
juncs = obI.nameProps.juncs;

%{


rGraph = zeros(length(sortPre2),length(sortPost2));
for i = 1:size(edges,1)
   targPre = find(sortPre2 == edges(i,2) );
   targPost = find(sortPost2 == edges(i,1));
   if ~isempty(targPre) & ~isempty(targPost);
      rGraph(targPre,targPost) = rGraph(targPre,targPost)+1; 
   end
end
sum(rGraph(:))

image(rGraph * 10)

%}

postJunc = zeros(length(sortPost2));
for i = 1:size(juncs,1)
   targPre = find(sortPost2 == juncs(i,2) );
   targPost = find(sortPost2 == juncs(i,1));
   if ~isempty(targPre) & ~isempty(targPost);
      postJunc(targPre,targPost) = postJunc(targPre,targPost)+1; 
      postJunc(targPost,targPre) = postJunc(targPost,targPre)+1; 

   else
       'missed'
   end
end
sum(postJunc(:))

image(postJunc * 10)


showCol(:,:,2) = postGraph2*10;
showCol(:,:,1) = postJunc *200;
showCol(:,:,3) = postGraph2*30;

image(uint8(showCol))


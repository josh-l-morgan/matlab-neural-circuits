

%Special use viewer (frustrated with fiji)
%{
Directions:
    1) Select directory with images
    2) Select command 
    3) Press return to apply command.  (can hold button for fast application)
Commands:
    d = forward
    a = backward
    w = brighter
    s = dimmer
    q = get next file
%}
while 1
clear all
TPN = GetMyDir;

'reading data'
TPNd = dir(TPN); TPNd = TPNd(3:length(TPNd));
TiffNames={};
for i = 1: size(TPNd,1)
    siz=length(TPNd(i).name);
    if TPNd(i).name(siz-2:siz)== 'tif'
        if TPNd(i).name(siz-8:siz-7)~='-R'
            TiffNames(length(TiffNames)+1,1)={TPNd(i).name};
        end
    end
end

for i = length(TiffNames):-1:1
    Iraw(:,:,:,i) = imread([TPN TiffNames{i}]);
end

%%
b = 1;
p = 1;
g = 'd';
Iraw(:,:,2,:) = Iraw(:,:,2,:) *2;
image(Iraw(:,:,:,p)*b),pause(.01)
while 1
inp = input('direction = ','s')

if ~strcmp(inp,'')
    g = inp;
end
if g == 'a'
    p = max(1,p - 1);
elseif g == 'd'
    p = min(size(Iraw,4),p+1);
elseif g == 's'
    b = b * .8;
elseif g == 'w'
    b = b * 1.2;
elseif g == 'q'
    break
end
    image(Iraw(:,:,:,p)*b),pause(.01)

% w = waitforbuttonpress
% inp = CurrentCharacter
end
end
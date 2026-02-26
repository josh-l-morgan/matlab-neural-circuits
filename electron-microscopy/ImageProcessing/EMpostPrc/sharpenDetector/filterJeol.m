colormap gray(256)
[TFN TPN] = GetMyFile;

I = imread([TPN TFN]);

%%
If = double(I);
%If = If * -1 + 255;
image(If)
subplot(2,2,1)
image(uint8(255-If(2000:2250,2000:2250)))
subplot(2,2,1)
%image(I)

% Sharpening Filter

%If = imfilter(If,[1;1;1]/3);

sFilt = [-1 -1 1 6 1 -1 -1 ];
sFilt = sFilt/4;
%sFilt = [ -1 3 -1 ] * 1;
If = imfilter(If,sFilt);
%If = imfilter(If,sFilt);
subplot(2,2,2)
image(uint8(255-If(2000:2250,2000:2250)))
%%


%If = imfilter(If,fspecial('disk',21));

% 
% 
If = imfilter(If,sFilt);
%If = imfilter(If,[1 ; 1; 1]/3);
%If = imfilter(If,[1 2 1 ; 2 10 2; 1 2 1]/25);
%If = medfilt2(If,[3 1]);
subplot(2,2,3)
image(uint8(255-If(2000:2250,2000:2250)))

%%
meanBright = mean(If,2);
for i = 1: size(If,1)
    Ivar(i) = var(If(i,:));
end
%bar(mean(I,1))
for i = 1:size(If,1)
    samp = Ivar(max(1,i-2):min(length(Ivar),i+2));
    meanSamp = mean(samp);
    If(i,:) = If(i,:) * meanSamp / Ivar(i);
end
subplot(2,2,4)
image(uint8(255-If(2000:2250,2000:2250)))

%%
imwrite(uint8(255-If),[TPN 'Filt' TFN],'tif','Compression','none')
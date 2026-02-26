clf
fsize = 200;
subplot(1,1,1)
xlim([1 fsize])
ylim([1 fsize])
for i = 1:100000000
%     spacer = repmat(' ',[1 ceil(rand*1000)]);
%     pause(.05)
%     disp(sprintf('%s%s',spacer,'stop interupting'))
    showText = 'stop interupting';
    textSize = length(showText)
    text(ceil(rand*fsize-textSize), ceil(rand*fsize-textSize),showText)
    pause(.1)
end
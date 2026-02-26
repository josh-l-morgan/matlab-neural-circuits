k=1;
[m n] = size(ONBTSLOW(k).spikes);
yMatrix = zeros(n, 1000);

for i=1:n
    input(i,:) = ONBTSLOW(k).spikes(i).input;
end

xBound = mean(input);
x = linspace(xBound(1), xBound(end), 1000);

for i=1:n
    fit = ONBTSLOW(k).spikes(i).rateFit;
    p = normcdf ( (fit(1) * x + fit(2)), 0, 1 );
    yMatrix(i,:) = p * fit(3) + fit (4);
end
y = mean(yMatrix);
yplusSEM = y + std(yMatrix)/sqrt(n);
yminusSEM = y - std(yMatrix)/sqrt(n);
plot(x,yplusSEM,'g', 'LineWidth', 1)
hold on
plot(x,yminusSEM,'g', 'LineWidth', 1)
plot(x,y,'g', 'LineWidth', 2)


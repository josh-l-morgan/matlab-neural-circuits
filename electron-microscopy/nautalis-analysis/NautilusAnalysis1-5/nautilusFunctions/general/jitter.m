function[j] = jitter(x,m)

%% and m evenly distributed random noise to x (m is width of noise not radius)

j = x + rand(size(x))* m - m/2;


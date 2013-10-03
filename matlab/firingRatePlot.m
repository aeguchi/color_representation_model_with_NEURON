tmp = zeros(30,30);
L5_000 = load('L5_FR_1000.txt');
for index=1:900
    x = mod(index-1,30)+1;
    y = ceil(index/30);
    tmp(y,x) = L5_000(index);
end
imagesc(tmp);
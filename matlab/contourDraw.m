FR = load('L5_FR_2000.txt');
FR_decomp = zeros(8,30,31);
figure;
number = 2;
set(gca, 'color', [0.5 0.5 0.5])
hold on

FR(find(FR>6))=6;

for col = 2:8
    FR_decomp(col,:,1:30) = reshape(FR(col,:),1,30,30);
    FR_decomp(col,1,31) = 6;
%     subplot(2,4,col);
    r = mod(ceil(col/4)+1,2);
    g = mod(ceil(col/2)+1,2);
    b = mod(ceil(col/1)+1,2);
%     if col==2
%         number = 1
%     else
%         number = 1
%     end
number = 1
    [C H] = contour(reshape(FR_decomp(col,:,:),30,31),number, 'LineColor',[r g b]);
    set (H, 'LineWidth', 2)%, 'LineStyle', ':'); 
    xlim([1 30]);
%     [C H] = contour(reshape(FR_decomp(col,:,:),30,30),1, 'LineColor',[r g b]);
%     set (H, 'LineWidth', 2); 
end




% sqd = zeros(300,3);
% div = size(img2,2)/300
% for i=1:300
%     sqd(i,:) = mean(img2(1,max([ceil(div*(i-1+0.1))-10 1]):min([ceil(div*i)+10 691]),:),2);
% end
% sqd2 = sqd./255

% sqd2 = load('wave400_700.txt')./255;
sqd2 = load('wave380_750.txt')./255;
randRGB = rand(10000,3);
%randRGB = load('colors_inputs.txt');

randWavelen = zeros(size(randRGB,1),1);
for i=1:size(randRGB,1)
    tmp = zeros(size(sqd2,1),3);
    tmp(:,1) = tmp(:,1)+randRGB(i,1);
    tmp(:,2) = tmp(:,2)+randRGB(i,2);
    tmp(:,3) = tmp(:,3)+randRGB(i,3);
    
    if((randRGB(i,1)+randRGB(i,2)+randRGB(i,3))>1.5)
        %tmp2 = abs(sqd2(20:350,1)-tmp(20:350,1))+abs(sqd2(20:350,2)-tmp(20:350,2))+abs(sqd2(20:350,3)-tmp(20:350,3));
        tmp2 = abs(sqd2(:,1)-tmp(:,1))+abs(sqd2(:,2)-tmp(:,2))+abs(sqd2(:,3)-tmp(:,3));
        %[find(tmp2==min(tmp2),1) min(tmp2)]
        randWavelen(i) = find(tmp2==min(tmp2),1,'last')+380;
    else
        randWavelen(i) = -1;
    end
end

figure();
subplot(9,2,1:2:15)
[counts, bins] =hist(randWavelen(find(randWavelen~=-1)),380:10:750);
counts = smooth(counts,0.3)
%bar(bins(2:size(bins,2)), counts(2:size(counts,2))/900)
plot(bins, counts/size(randRGB,1),'LineWidth', 2, 'color',[0.8 0.8 0.8], 'LineStyle', '--')
hold on

subplot(9,2,2:2:16)
plot(bins, counts/size(randRGB,1),'LineWidth', 2, 'color',[0.8 0.8 0.8], 'LineStyle', '--')
hold on




% sqd2 = load('wave400_700.txt')./255;
sqd2 = load('wave380_750.txt')./255;
randRGB = sqrt(load('colors_half_u.txt'));
%randRGB = load('colors_inputs.txt');

randWavelen = zeros(size(randRGB,1),1);

for i=1:size(randRGB,1)
    tmp = zeros(size(sqd2,1),3);
    tmp(:,1) = tmp(:,1)+randRGB(i,1);
    tmp(:,2) = tmp(:,2)+randRGB(i,2);
    tmp(:,3) = tmp(:,3)+randRGB(i,3);
    
    if((randRGB(i,1)+randRGB(i,2)+randRGB(i,3))>1.5)
        %tmp2 = abs(sqd2(20:350,1)-tmp(20:350,1))+abs(sqd2(20:350,2)-tmp(20:350,2))+abs(sqd2(20:350,3)-tmp(20:350,3));
        tmp2 = abs(sqd2(:,1)-tmp(:,1))+abs(sqd2(:,2)-tmp(:,2))+abs(sqd2(:,3)-tmp(:,3));
        %[find(tmp2==min(tmp2),1) min(tmp2)]
        randWavelen(i) = find(tmp2==min(tmp2),1,'last')+380;
    else
        randWavelen(i) = -1;
    end
end

subplot(9,2,1:2:16)
[counts, bins] =hist(randWavelen(find(randWavelen~=-1)),380:10:750);
%bar(bins(2:size(bins,2)), counts(2:size(counts,2))/900)
plot(bins, counts/size(randRGB,1),'LineWidth', 2, 'color','k')
xlim([380 750])
ylabel({'Frequency of cells';'sensitive to different wavelength [freq/10nm]'})
%ylim([0 0.05])
set(gca,'xtick',[])
subplot(9,2,17)
image(reshape(sqd2,1,370,3))
set(gca,'ytick',[])
% set(gca,'XTickLabel',(get(gca,'XTickLabel')+400));
xlim([0 370])
xlabel('Wavelength [nm]')
set(gca,'xticklabel',get(gca,'xtick')+380)





%randRGB = rand(900,3);
%randRGB = sqrt(load('colors_t.txt'));
randRGB = sqrt(load('colors_half_t.txt'));
%randRGB = sqrt(load('colors_nat_t.txt'));


randWavelen = zeros(size(randRGB,1),1);

for i=1:size(randRGB,1)
    tmp = zeros(size(sqd2,1),3);
    tmp(:,1) = tmp(:,1)+randRGB(i,1);
    tmp(:,2) = tmp(:,2)+randRGB(i,2);
    tmp(:,3) = tmp(:,3)+randRGB(i,3);
    
    if((randRGB(i,1)+randRGB(i,2)+randRGB(i,3))>1.5)
        tmp2 = abs(sqd2(:,1)-tmp(:,1))+abs(sqd2(:,2)-tmp(:,2))+abs(sqd2(:,3)-tmp(:,3));
        %tmp2 = abs(sqd2(:,1)-tmp(:,1))+abs(sqd2(:,2)-tmp(:,2))+abs(sqd2(:,3)-tmp(:,3));
        %[find(tmp2==min(tmp2),1) min(tmp2)]
        randWavelen(i) = find(tmp2==min(tmp2),1,'last')+350;
    else
        randWavelen(i) = -1;
    end
end

subplot(9,2,2:2:16)
[counts, bins] =hist(randWavelen(find(randWavelen~=-1)),380:10:750);
%bar(bins(2:size(bins,2)), counts(2:size(counts,2))/900)
plot(bins, counts/900,'LineWidth', 2, 'color','k')
xlim([380 750])
ylim([0 0.05])
ylabel({'Frequency of cells';'sensitive to different wavelength [freq/10nm]'})
set(gca,'xtick',[])
subplot(9,2,18)
image(reshape(sqd2,1,370,3))
xlim([0 370])
set(gca,'ytick',[])
xlabel('Wavelength [nm]')
set(gca,'xticklabel',get(gca,'xtick')+380)



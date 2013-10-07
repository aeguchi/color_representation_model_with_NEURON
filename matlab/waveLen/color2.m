

% sqd = zeros(370,3);
% div = size(avg,2)/370
% for i=1:370
%     sqd(i,:) = mean(avg(1,ceil(div*(i-1+0.1)):ceil(div*i),:),2)
% end



sqd2 = load('wave380_750.txt')./255;
%randRGB = rand(900,3);
%randRGB = load('colors_u.txt');
randRGB = load('colors_half_u.txt');

figure
% subplot(3,2,1)
[counts, bins] =hist(sqrt(randRGB(:,1)),0:0.05:1);
%bar(bins, counts/900)
%h = findobj(gca,'Type','patch');
%set(h,'FaceColor','r')

subplot(1,2,1)
p = plot(bins, counts/900,'color','r', 'LineWidth', 2)
hold on
[counts, bins] =hist(sqrt(randRGB(:,2)),0:0.05:1);
plot(bins, counts/900,'color','g', 'LineWidth', 2, 'LineStyle', '--')
[counts, bins] =hist(sqrt(randRGB(:,3)),0:0.05:1);
plot(bins, counts/900,'color','b', 'LineWidth', 2, 'LineStyle', ':')
xlim([0 1])
ylim([0 0.35])
xlabel('sensitivity')
ylabel({'Frequency of cells';'sensitive to each RGB [freq/0.05]'})
% ylabel({'number of neurons';'sensitive to Red'})


% subplot(3,2,3)
% [counts, bins] =hist(randRGB(:,2),0:0.05:1);
% bar(bins, counts/900)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','g')
% xlim([0 1])
% ylim([0 0.3])
% xlabel('sensitivity')
% ylabel({'number of neurons';'sensitive to Green'})
% 
% subplot(3,2,5)
% [counts, bins] =hist(randRGB(:,3),0:0.05:1);
% bar(bins, counts/900)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','b')
% xlim([0 1])
% ylim([0 0.3])
% xlabel('sensitivity')
% ylabel({'number of neurons';'sensitive to Blue'})



%randRGB = load('colors_t.txt');
randRGB = load('colors_half_t.txt');

% subplot(3,2,2)
[counts, bins] =hist(sqrt(randRGB(:,1)),0:0.05:1);
% bar(bins, counts/900)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r')
subplot(1,2,2)
plot(bins, counts/900,'color','r', 'LineWidth', 2)
hold on
[counts, bins] =hist(sqrt(randRGB(:,2)),0:0.05:1);
plot(bins, counts/900,'color','g', 'LineWidth', 2, 'LineStyle', '--')
[counts, bins] =hist(sqrt(randRGB(:,3)),0:0.05:1);
plot(bins, counts/900,'color','b', 'LineWidth', 2, 'LineStyle', ':')
xlim([0 1])
ylim([0 0.5])
xlabel('sensitivity')
ylabel({'number of neurons';'sensitive to Red'})


xlim([0 1])
ylim([0 0.35])
xlabel('sensitivity')
ylabel({'Frequency of cells';'sensitive to each RGB [freq/0.05]'})
% ylabel({'number of neurons';'sensitive to Red'})

% subplot(3,2,4)
% [counts, bins] =hist(randRGB(:,2),0:0.05:1);
% bar(bins, counts/900)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','g')
% xlim([0 1])
% ylim([0 0.3])
% xlabel('sensitivity')
% ylabel({'number of neurons';'sensitive to Green'})
% 
% subplot(3,2,6)
% [counts, bins] =hist(randRGB(:,3),0:0.05:1);
% bar(bins, counts/900)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','b')
% xlim([0 1])
% ylim([0 0.3])
% xlabel('sensitivity')
% ylabel({'number of neurons';'sensitive to Blue'})

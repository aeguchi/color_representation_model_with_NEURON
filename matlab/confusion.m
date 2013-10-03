img(1,:,:) = [[0 0 0];[0 0 1];[0 1 0];[0 1 1];[1 0 0];[1 0 1];[1 1 0];[1 1 1]];
img2 = reshape([[0 0 0];[0 0 1];[0 1 0];[0 1 1];[1 0 0];[1 0 1];[1 1 0];[1 1 1]],8,1,3);

subplot(10,10,1:10:81)
image(img2)
set(gca,'YTick',[])
set(gca,'XTick',[])
subplot(10,10,92:100)
image(img)
set(gca,'XTick',[])
set(gca,'YTick',[])

subplot(10,10,[2:10 12:20 22:30 32:40 42:50 52:60 62:70 72:80 82:90])
imagesc(untrained,[0 1/8])
set(gca,'XTick',[])
set(gca,'YTick',[])

imagesc(trained,[0 1/8])
set(gca,'YTick',[])
set(gca,'XTick',[])
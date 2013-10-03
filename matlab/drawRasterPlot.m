nDim = 30
tstop = 300
dt = 0.025
img = zeros(1, 40, 3);
img(1,:,:) = [[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0.5 0.5 0.5];
              [0 0 1];[0 0 1];[0 0 1];[0 0 1];[0.5 0.5 0.5];
              [0 1 0];[0 1 0];[0 1 0];[0 1 0];[0.5 0.5 0.5];
              [0 1 1];[0 1 1];[0 1 1];[0 1 1];[0.5 0.5 0.5];
              [1 0 0];[1 0 0];[1 0 0];[1 0 0];[0.5 0.5 0.5];
              [1 0 1];[1 0 1];[1 0 1];[1 0 1];[0.5 0.5 0.5];
              [1 1 0];[1 1 0];[1 1 0];[1 1 0];[0.5 0.5 0.5];
              [1 1 1];[1 1 1];[1 1 1];[1 1 1];[0.5 0.5 0.5];];



L23000 = load('L23_firing_0_[0, 0, 0].txt');
L23000_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23000_sq(:,t+1) = mean(L23000(:,index:index+40),2);
end
clear L23000
L23001 = load('L23_firing_0_[0, 0, 1].txt');
L23001_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23001_sq(:,t+1) = mean(L23001(:,index:index+40),2);
end
clear L23001
L23010 = load('L23_firing_0_[0, 1, 0].txt');
L23010_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23010_sq(:,t+1) = mean(L23010(:,index:index+40),2);
end
clear L23010
L23011 = load('L23_firing_0_[0, 1, 1].txt');
L23011_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23011_sq(:,t+1) = mean(L23011(:,index:index+40),2);
end
clear L23011
L23100 = load('L23_firing_0_[1, 0, 0].txt');
L23100_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23100_sq(:,t+1) = mean(L23100(:,index:index+40),2);
end
clear L23100
L23101 = load('L23_firing_0_[1, 0, 1].txt');
L23101_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23101_sq(:,t+1) = mean(L23101(:,index:index+40),2);
end
clear L23101
L23110 = load('L23_firing_0_[1, 1, 0].txt');
L23110_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23110_sq(:,t+1) = mean(L23110(:,index:index+40),2);
end
clear L23110
L23111 = load('L23_firing_0_[1, 1, 1].txt');
L23111_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23111_sq(:,t+1) = mean(L23111(:,index:index+40),2);
end
clear L23111
comb_L23_0 = [L23000_sq L23001_sq L23010_sq L23011_sq L23100_sq L23101_sq L23110_sq L23111_sq];




L4000 = load('L4_firing_0_[0, 0, 0].txt');
L4000_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4000_sq(:,t+1) = mean(L4000(:,index:index+40),2);
end
clear L4000
L4001 = load('L4_firing_0_[0, 0, 1].txt');
L4001_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4001_sq(:,t+1) = mean(L4001(:,index:index+40),2);
end
clear L4001
L4010 = load('L4_firing_0_[0, 1, 0].txt');
L4010_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4010_sq(:,t+1) = mean(L4010(:,index:index+40),2);
end
clear L4010
L4011 = load('L4_firing_0_[0, 1, 1].txt');
L4011_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4011_sq(:,t+1) = mean(L4011(:,index:index+40),2);
end
clear L4011
L4100 = load('L4_firing_0_[1, 0, 0].txt');
L4100_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4100_sq(:,t+1) = mean(L4100(:,index:index+40),2);
end
clear L4100
L4101 = load('L4_firing_0_[1, 0, 1].txt');
L4101_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4101_sq(:,t+1) = mean(L4101(:,index:index+40),2);
end
clear L4101
L4110 = load('L4_firing_0_[1, 1, 0].txt');
L4110_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4110_sq(:,t+1) = mean(L4110(:,index:index+40),2);
end
clear L4110
L4111 = load('L4_firing_0_[1, 1, 1].txt');
L4111_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4111_sq(:,t+1) = mean(L4111(:,index:index+40),2);
end
clear L4111
comb_L4_0 = [L4000_sq L4001_sq L4010_sq L4011_sq L4100_sq L4101_sq L4110_sq L4111_sq];






L5000 = load('L5_firing_0_[0, 0, 0].txt');
L5000_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5000_sq(:,t+1) = mean(L5000(:,index:index+40),2);
end
clear L5000
L5001 = load('L5_firing_0_[0, 0, 1].txt');
L5001_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5001_sq(:,t+1) = mean(L5001(:,index:index+40),2);
end
clear L5001
L5010 = load('L5_firing_0_[0, 1, 0].txt');
L5010_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5010_sq(:,t+1) = mean(L5010(:,index:index+40),2);
end
clear L5010
L5011 = load('L5_firing_0_[0, 1, 1].txt');
L5011_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5011_sq(:,t+1) = mean(L5011(:,index:index+40),2);
end
clear L5011
L5100 = load('L5_firing_0_[1, 0, 0].txt');
L5100_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5100_sq(:,t+1) = mean(L5100(:,index:index+40),2);
end
clear L5100
L5101 = load('L5_firing_0_[1, 0, 1].txt');
L5101_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5101_sq(:,t+1) = mean(L5101(:,index:index+40),2);
end
clear L5101
L5110 = load('L5_firing_0_[1, 1, 0].txt');
L5110_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5110_sq(:,t+1) = mean(L5110(:,index:index+40),2);
end
clear L5110
L5111 = load('L5_firing_0_[1, 1, 1].txt');
L5111_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5111_sq(:,t+1) = mean(L5111(:,index:index+40),2);
end
clear L5111
comb_L5_0 = [L5000_sq L5001_sq L5010_sq L5011_sq L5100_sq L5101_sq L5110_sq L5111_sq];




%%%input channel load%%
itr = '0'
L000_sq = zeros(nDim*nDim,tstop+1);
L001_sq = zeros(nDim*nDim,tstop+1);
L010_sq = zeros(nDim*nDim,tstop+1);
L011_sq = zeros(nDim*nDim,tstop+1);
L100_sq = zeros(nDim*nDim,tstop+1);
L101_sq = zeros(nDim*nDim,tstop+1);
L110_sq = zeros(nDim*nDim,tstop+1);
L111_sq = zeros(nDim*nDim,tstop+1);
L000 = load(['L_firing_' itr '_[0, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L000_sq(:,t+1) = mean(L000(:,index:index+40),2);
end
clear L000;
L001 = load(['L_firing_' itr '_[0, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L001_sq(:,t+1) = mean(L001(:,index:index+40),2);
end
clear L001;
L010 = load(['L_firing_' itr '_[0, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L010_sq(:,t+1) = mean(L010(:,index:index+40),2);
end
clear L010;
L011 = load(['L_firing_' itr '_[0, 1, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L011_sq(:,t+1) = mean(L011(:,index:index+40),2);
end
clear L011;
L100 = load(['L_firing_' itr '_[1, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L100_sq(:,t+1) = mean(L100(:,index:index+40),2);
end
clear L100;
L101 = load(['L_firing_' itr '_[1, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L101_sq(:,t+1) = mean(L101(:,index:index+40),2);
end
clear L101;
L110 = load(['L_firing_' itr '_[1, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L110_sq(:,t+1) = mean(L110(:,index:index+40),2);
end
clear L110;
L111 = load(['L_firing_' itr '_[1, 1, 1].txt']);

for t = 0:tstop-1
    index = 2+(t)*40;
    L111_sq(:,t+1) = mean(L111(:,index:index+40),2);
end
clear L111;
comb_L = [L000_sq L001_sq L010_sq L011_sq L100_sq L101_sq L110_sq L111_sq];


C1000_sq = zeros(nDim*nDim,tstop+1);
C1001_sq = zeros(nDim*nDim,tstop+1);
C1010_sq = zeros(nDim*nDim,tstop+1);
C1011_sq = zeros(nDim*nDim,tstop+1);
C1100_sq = zeros(nDim*nDim,tstop+1);
C1101_sq = zeros(nDim*nDim,tstop+1);
C1110_sq = zeros(nDim*nDim,tstop+1);
C1111_sq = zeros(nDim*nDim,tstop+1);
C1000 = load(['C1_firing_' itr '_[0, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C1000_sq(:,t+1) = mean(C1000(:,index:index+40),2);
end
clear C1000;
C1001 = load(['C1_firing_' itr '_[0, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C1001_sq(:,t+1) = mean(C1001(:,index:index+40),2);
end
clear C1001;
C1010 = load(['C1_firing_' itr '_[0, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C1010_sq(:,t+1) = mean(C1010(:,index:index+40),2);
end
clear C1010;
C1011 = load(['C1_firing_' itr '_[0, 1, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C1011_sq(:,t+1) = mean(C1011(:,index:index+40),2);
end
clear C1011;
C1100 = load(['C1_firing_' itr '_[1, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C1100_sq(:,t+1) = mean(C1100(:,index:index+40),2);
end
clear C1100;
C1101 = load(['C1_firing_' itr '_[1, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C1101_sq(:,t+1) = mean(C1101(:,index:index+40),2);
end
clear C1101;
C1110 = load(['C1_firing_' itr '_[1, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C1110_sq(:,t+1) = mean(C1110(:,index:index+40),2);
end
clear C1110;
C1111 = load(['C1_firing_' itr '_[1, 1, 1].txt']);

for t = 0:tstop-1
    index = 2+(t)*40;
    C1111_sq(:,t+1) = mean(C1111(:,index:index+40),2);
end
clear C1111;
comb_C1 = [C1000_sq C1001_sq C1010_sq C1011_sq C1100_sq C1101_sq C1110_sq C1111_sq];


C2000_sq = zeros(nDim*nDim,tstop+1);
C2001_sq = zeros(nDim*nDim,tstop+1);
C2010_sq = zeros(nDim*nDim,tstop+1);
C2011_sq = zeros(nDim*nDim,tstop+1);
C2100_sq = zeros(nDim*nDim,tstop+1);
C2101_sq = zeros(nDim*nDim,tstop+1);
C2110_sq = zeros(nDim*nDim,tstop+1);
C2111_sq = zeros(nDim*nDim,tstop+1);
C2000 = load(['C2_firing_' itr '_[0, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C2000_sq(:,t+1) = mean(C2000(:,index:index+40),2);
end
clear C2000;
C2001 = load(['C2_firing_' itr '_[0, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C2001_sq(:,t+1) = mean(C2001(:,index:index+40),2);
end
clear C2001;
C2010 = load(['C2_firing_' itr '_[0, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C2010_sq(:,t+1) = mean(C2010(:,index:index+40),2);
end
clear C2010;
C2011 = load(['C2_firing_' itr '_[0, 1, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C2011_sq(:,t+1) = mean(C2011(:,index:index+40),2);
end
clear C2011;
C2100 = load(['C2_firing_' itr '_[1, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C2100_sq(:,t+1) = mean(C2100(:,index:index+40),2);
end
clear C2100;
C2101 = load(['C2_firing_' itr '_[1, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C2101_sq(:,t+1) = mean(C2101(:,index:index+40),2);
end
clear C2101;
C2110 = load(['C2_firing_' itr '_[1, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    C2110_sq(:,t+1) = mean(C2110(:,index:index+40),2);
end
clear C2110;
C2111 = load(['C2_firing_' itr '_[1, 1, 1].txt']);

for t = 0:tstop-1
    index = 2+(t)*40;
    C2111_sq(:,t+1) = mean(C2111(:,index:index+40),2);
end
clear C2111;
comb_C2 = [C2000_sq C2001_sq C2010_sq C2011_sq C2100_sq C2101_sq C2110_sq C2111_sq];


%combined 
fig = figure
set(fig, 'Position', [0 0 560 720])
subplot(19,1,1:6)
Input = comb_L;
imagesc((1-Input))
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylabel('L channel (R+G) (900 neurons)')

subplot(19,1,7:12)
Input = comb_C1;
imagesc((1-Input))
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylabel('C1 channel (R-G) (900 neurons)')

subplot(19,1,13:18)
Input = comb_C2;
imagesc((1-Input))
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylabel('C2 channel (R+G-B) (900 neurons)')


subplot(19,1,19)
image(img);
set(gca, 'XTick', [])
set(gca, 'YTick', [])
xlabel('input colors [240 ms/color] with 60ms intervals');













%000
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23000_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)
%saveas(fig, 'r_000_u.png');

%001
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23001_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)

%010
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23010_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)

%011
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23011_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)

%100
Input = L23100_sq;
fig = figure
set(fig, 'Position', [0 0 560 720])
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)

%101
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23101_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)


%110
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23110_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)


%111
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23111_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)


%combined 
fig = figure
set(fig, 'Position', [0 0 560 720])
subplot(19,1,1:5)
Input = [comb_L5_0;comb_L23_0;comb_L4_0];
imagesc((1-Input(1:nDim*nDim,:)))
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylabel('V1_5 (900 cells)');

subplot(19,1,6)
bar(sum(Input(1:nDim*nDim,:),1))
xlim([0 (tstop+1)*8-1])
set(gca, 'XTick', [])

subplot(19,1,7:11)
imagesc((1-Input(nDim*nDim+1:nDim*nDim*2,:)))
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylabel('V1_2_/_3 (900 cells)');

subplot(19,1,12)
bar(sum(Input(nDim*nDim+1:nDim*nDim*2,:),1))
xlim([0 (tstop+1)*8-1])
set(gca, 'XTick', [])

subplot(19,1,13:17)
imagesc((1-Input(nDim*nDim*2+1:nDim*nDim*3,:)))
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylabel('V1_4 (900 cells)');

subplot(19,1,18)
bar(sum(Input(nDim*nDim*2+1:nDim*nDim*3,:),1))
xlim([0 (tstop+1)*8-1])
set(gca, 'XTick', [])

subplot(19,1,19)
image(img);
set(gca, 'XTick', [])
set(gca, 'YTick', [])
xlabel('input colors [240 ms/color] with 60ms intervals');





























itr = '2000';

L23000_2_sq = zeros(nDim*nDim,tstop+1);
L23001_2_sq = zeros(nDim*nDim,tstop+1);
L23010_2_sq = zeros(nDim*nDim,tstop+1);
L23011_2_sq = zeros(nDim*nDim,tstop+1);
L23100_2_sq = zeros(nDim*nDim,tstop+1);
L23101_2_sq = zeros(nDim*nDim,tstop+1);
L23110_2_sq = zeros(nDim*nDim,tstop+1);
L23111_2_sq = zeros(nDim*nDim,tstop+1);
L23000_2 = load(['L23_firing_' itr '_[0, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23000_2_sq(:,t+1) = mean(L23000_2(:,index:index+40),2);
end
clear L23000_2;
L23001_2 = load(['L23_firing_' itr '_[0, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23001_2_sq(:,t+1) = mean(L23001_2(:,index:index+40),2);
end
clear L23001_2;
L23010_2 = load(['L23_firing_' itr '_[0, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23010_2_sq(:,t+1) = mean(L23010_2(:,index:index+40),2);
end
clear L23010_2;
L23011_2 = load(['L23_firing_' itr '_[0, 1, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23011_2_sq(:,t+1) = mean(L23011_2(:,index:index+40),2);
end
clear L23011_2;
L23100_2 = load(['L23_firing_' itr '_[1, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23100_2_sq(:,t+1) = mean(L23100_2(:,index:index+40),2);
end
clear L23100_2;
L23101_2 = load(['L23_firing_' itr '_[1, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23101_2_sq(:,t+1) = mean(L23101_2(:,index:index+40),2);
end
clear L23101_2;
L23110_2 = load(['L23_firing_' itr '_[1, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23110_2_sq(:,t+1) = mean(L23110_2(:,index:index+40),2);
end
clear L23110_2;
L23111_2 = load(['L23_firing_' itr '_[1, 1, 1].txt']);

for t = 0:tstop-1
    index = 2+(t)*40;
    L23111_2_sq(:,t+1) = mean(L23111_2(:,index:index+40),2);
end
clear L23111_2;
comb_L23_2 = [L23000_2_sq L23001_2_sq L23010_2_sq L23011_2_sq L23100_2_sq L23101_2_sq L23110_2_sq L23111_2_sq];



L4000_2_sq = zeros(nDim*nDim,tstop+1);
L4001_2_sq = zeros(nDim*nDim,tstop+1);
L4010_2_sq = zeros(nDim*nDim,tstop+1);
L4011_2_sq = zeros(nDim*nDim,tstop+1);
L4100_2_sq = zeros(nDim*nDim,tstop+1);
L4101_2_sq = zeros(nDim*nDim,tstop+1);
L4110_2_sq = zeros(nDim*nDim,tstop+1);
L4111_2_sq = zeros(nDim*nDim,tstop+1);
L4000_2 = load(['L4_firing_' itr '_[0, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4000_2_sq(:,t+1) = mean(L4000_2(:,index:index+40),2);
end
clear L4000_2;
L4001_2 = load(['L4_firing_' itr '_[0, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4001_2_sq(:,t+1) = mean(L4001_2(:,index:index+40),2);
end
clear L4001_2;
L4010_2 = load(['L4_firing_' itr '_[0, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4010_2_sq(:,t+1) = mean(L4010_2(:,index:index+40),2);
end
clear L4010_2;
L4011_2 = load(['L4_firing_' itr '_[0, 1, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4011_2_sq(:,t+1) = mean(L4011_2(:,index:index+40),2);
end
clear L4011_2;
L4100_2 = load(['L4_firing_' itr '_[1, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4100_2_sq(:,t+1) = mean(L4100_2(:,index:index+40),2);
end
clear L4100_2;
L4101_2 = load(['L4_firing_' itr '_[1, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4101_2_sq(:,t+1) = mean(L4101_2(:,index:index+40),2);
end
clear L4101_2;
L4110_2 = load(['L4_firing_' itr '_[1, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L4110_2_sq(:,t+1) = mean(L4110_2(:,index:index+40),2);
end
clear L4110_2;
L4111_2 = load(['L4_firing_' itr '_[1, 1, 1].txt']);

for t = 0:tstop-1
    index = 2+(t)*40;
    L4111_2_sq(:,t+1) = mean(L4111_2(:,index:index+40),2);
end
clear L4111_2;
comb_L4_2 = [L4000_2_sq L4001_2_sq L4010_2_sq L4011_2_sq L4100_2_sq L4101_2_sq L4110_2_sq L4111_2_sq];


L5000_2_sq = zeros(nDim*nDim,tstop+1);
L5001_2_sq = zeros(nDim*nDim,tstop+1);
L5010_2_sq = zeros(nDim*nDim,tstop+1);
L5011_2_sq = zeros(nDim*nDim,tstop+1);
L5100_2_sq = zeros(nDim*nDim,tstop+1);
L5101_2_sq = zeros(nDim*nDim,tstop+1);
L5110_2_sq = zeros(nDim*nDim,tstop+1);
L5111_2_sq = zeros(nDim*nDim,tstop+1);
L5000_2 = load(['L5_firing_' itr '_[0, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5000_2_sq(:,t+1) = mean(L5000_2(:,index:index+40),2);
end
clear L5000_2;
L5001_2 = load(['L5_firing_' itr '_[0, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5001_2_sq(:,t+1) = mean(L5001_2(:,index:index+40),2);
end
clear L5001_2;
L5010_2 = load(['L5_firing_' itr '_[0, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5010_2_sq(:,t+1) = mean(L5010_2(:,index:index+40),2);
end
clear L5010_2;
L5011_2 = load(['L5_firing_' itr '_[0, 1, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5011_2_sq(:,t+1) = mean(L5011_2(:,index:index+40),2);
end
clear L5011_2;
L5100_2 = load(['L5_firing_' itr '_[1, 0, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5100_2_sq(:,t+1) = mean(L5100_2(:,index:index+40),2);
end
clear L5100_2;
L5101_2 = load(['L5_firing_' itr '_[1, 0, 1].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5101_2_sq(:,t+1) = mean(L5101_2(:,index:index+40),2);
end
clear L5101_2;
L5110_2 = load(['L5_firing_' itr '_[1, 1, 0].txt']);
for t = 0:tstop-1
    index = 2+(t)*40;
    L5110_2_sq(:,t+1) = mean(L5110_2(:,index:index+40),2);
end
clear L5110_2;
L5111_2 = load(['L5_firing_' itr '_[1, 1, 1].txt']);

for t = 0:tstop-1
    index = 2+(t)*40;
    L5111_2_sq(:,t+1) = mean(L5111_2(:,index:index+40),2);
end
clear L5111_2;
comb_L5_2 = [L5000_2_sq L5001_2_sq L5010_2_sq L5011_2_sq L5100_2_sq L5101_2_sq L5110_2_sq L5111_2_sq];









% for t = 0:tstop-1
%     index = 2+(t)*40;
%     L4000_2_sq(:,t+1) = mean(L4000_2(:,index:index+40),2);
%     L4001_2_sq(:,t+1) = mean(L4001_2(:,index:index+40),2);
%     L4010_2_sq(:,t+1) = mean(L4010_2(:,index:index+40),2);
%     L4011_2_sq(:,t+1) = mean(L4011_2(:,index:index+40),2);
%     L4100_2_sq(:,t+1) = mean(L4100_2(:,index:index+40),2);
%     L4101_2_sq(:,t+1) = mean(L4101_2(:,index:index+40),2);
%     L4110_2_sq(:,t+1) = mean(L4110_2(:,index:index+40),2);
%     L4111_2_sq(:,t+1) = mean(L4111_2(:,index:index+40),2);
% end
% L4000_2 = 0;
% L4001_2 = 0;
% L4010_2 = 0;
% L4011_2 = 0;
% L4100_2 = 0;
% L4101_2 = 0;
% L4110_2 = 0;
% L4111_2 = 0;
% comb_L4_2 = [L4000_2_sq L4001_2_sq L4010_2_sq L4011_2_sq L4100_2_sq L4101_2_sq L4110_2_sq L4111_2_sq];



% L000_2 = load(['L_firing_' itr '_[0, 0, 0].txt']);
% L001_2 = load(['L_firing_' itr '_[0, 0, 1].txt']);
% L010_2 = load(['L_firing_' itr '_[0, 1, 0].txt']);
% L011_2 = load(['L_firing_' itr '_[0, 1, 1].txt']);
% L100_2 = load(['L_firing_' itr '_[1, 0, 0].txt']);
% L101_2 = load(['L_firing_' itr '_[1, 0, 1].txt']);
% L110_2 = load(['L_firing_' itr '_[1, 1, 0].txt']);
% L111_2 = load(['L_firing_' itr '_[1, 1, 1].txt']);
% 
% L000_2_sq = zeros(nDim*nDim,tstop+1);
% L001_2_sq = zeros(nDim*nDim,tstop+1);
% L010_2_sq = zeros(nDim*nDim,tstop+1);
% L011_2_sq = zeros(nDim*nDim,tstop+1);
% L100_2_sq = zeros(nDim*nDim,tstop+1);
% L101_2_sq = zeros(nDim*nDim,tstop+1);
% L110_2_sq = zeros(nDim*nDim,tstop+1);
% L111_2_sq = zeros(nDim*nDim,tstop+1);
% for t = 0:tstop-1
%     index = 2+(t)*40;
%     L000_2_sq(:,t+1) = mean(L000_2(:,index:index+40),2);
%     L001_2_sq(:,t+1) = mean(L001_2(:,index:index+40),2);
%     L010_2_sq(:,t+1) = mean(L010_2(:,index:index+40),2);
%     L011_2_sq(:,t+1) = mean(L011_2(:,index:index+40),2);
%     L100_2_sq(:,t+1) = mean(L100_2(:,index:index+40),2);
%     L101_2_sq(:,t+1) = mean(L101_2(:,index:index+40),2);
%     L110_2_sq(:,t+1) = mean(L110_2(:,index:index+40),2);
%     L111_2_sq(:,t+1) = mean(L111_2(:,index:index+40),2);
% end
% clear L000_2 L001_2 L010_2 L011_2 L100_2 L101_2 L110_2 L111_2
% comb_L_2 = [L000_2_sq L001_2_sq L010_2_sq L011_2_sq L100_2_sq L101_2_sq L110_2_sq L111_2_sq];
% 
% 
% C1000_2 = load(['C1_firing_' itr '_[0, 0, 0].txt']);
% C1001_2 = load(['C1_firing_' itr '_[0, 0, 1].txt']);
% C1010_2 = load(['C1_firing_' itr '_[0, 1, 0].txt']);
% C1011_2 = load(['C1_firing_' itr '_[0, 1, 1].txt']);
% C1100_2 = load(['C1_firing_' itr '_[1, 0, 0].txt']);
% C1101_2 = load(['C1_firing_' itr '_[1, 0, 1].txt']);
% C1110_2 = load(['C1_firing_' itr '_[1, 1, 0].txt']);
% C1111_2 = load(['C1_firing_' itr '_[1, 1, 1].txt']);
% 
% C1000_2_sq = zeros(nDim*nDim,tstop+1);
% C1001_2_sq = zeros(nDim*nDim,tstop+1);
% C1010_2_sq = zeros(nDim*nDim,tstop+1);
% C1011_2_sq = zeros(nDim*nDim,tstop+1);
% C1100_2_sq = zeros(nDim*nDim,tstop+1);
% C1101_2_sq = zeros(nDim*nDim,tstop+1);
% C1110_2_sq = zeros(nDim*nDim,tstop+1);
% C1111_2_sq = zeros(nDim*nDim,tstop+1);
% for t = 0:tstop-1
%     index = 2+(t)*40;
%     C1000_2_sq(:,t+1) = mean(C1000_2(:,index:index+40),2);
%     C1001_2_sq(:,t+1) = mean(C1001_2(:,index:index+40),2);
%     C1010_2_sq(:,t+1) = mean(C1010_2(:,index:index+40),2);
%     C1011_2_sq(:,t+1) = mean(C1011_2(:,index:index+40),2);
%     C1100_2_sq(:,t+1) = mean(C1100_2(:,index:index+40),2);
%     C1101_2_sq(:,t+1) = mean(C1101_2(:,index:index+40),2);
%     C1110_2_sq(:,t+1) = mean(C1110_2(:,index:index+40),2);
%     C1111_2_sq(:,t+1) = mean(C1111_2(:,index:index+40),2);
% end
% clear C1000_2 C1001_2 C1010_2 C1011_2 C1100_2 C1101_2 C1110_2 C1111_2
% comb_C1_2 = [C1000_2_sq C1001_2_sq C1010_2_sq C1011_2_sq C1100_2_sq C1101_2_sq C1110_2_sq C1111_2_sq];
% 
% C2000_2 = load(['C2_firing_' itr '_[0, 0, 0].txt']);
% C2001_2 = load(['C2_firing_' itr '_[0, 0, 1].txt']);
% C2010_2 = load(['C2_firing_' itr '_[0, 1, 0].txt']);
% C2011_2 = load(['C2_firing_' itr '_[0, 1, 1].txt']);
% C2100_2 = load(['C2_firing_' itr '_[1, 0, 0].txt']);
% C2101_2 = load(['C2_firing_' itr '_[1, 0, 1].txt']);
% C2110_2 = load(['C2_firing_' itr '_[1, 1, 0].txt']);
% C2111_2 = load(['C2_firing_' itr '_[1, 1, 1].txt']);
% 
% C2000_2_sq = zeros(nDim*nDim,tstop+1);
% C2001_2_sq = zeros(nDim*nDim,tstop+1);
% C2010_2_sq = zeros(nDim*nDim,tstop+1);
% C2011_2_sq = zeros(nDim*nDim,tstop+1);
% C2100_2_sq = zeros(nDim*nDim,tstop+1);
% C2101_2_sq = zeros(nDim*nDim,tstop+1);
% C2110_2_sq = zeros(nDim*nDim,tstop+1);
% C2111_2_sq = zeros(nDim*nDim,tstop+1);
% for t = 0:tstop-1
%     index = 2+(t)*40;
%     C2000_2_sq(:,t+1) = mean(C2000_2(:,index:index+40),2);
%     C2001_2_sq(:,t+1) = mean(C2001_2(:,index:index+40),2);
%     C2010_2_sq(:,t+1) = mean(C2010_2(:,index:index+40),2);
%     C2011_2_sq(:,t+1) = mean(C2011_2(:,index:index+40),2);
%     C2100_2_sq(:,t+1) = mean(C2100_2(:,index:index+40),2);
%     C2101_2_sq(:,t+1) = mean(C2101_2(:,index:index+40),2);
%     C2110_2_sq(:,t+1) = mean(C2110_2(:,index:index+40),2);
%     C2111_2_sq(:,t+1) = mean(C2111_2(:,index:index+40),2);
% end
% clear C2000_2 C2001_2 C2010_2 C2011_2 C2100_2 C2101_2 C2110_2 C2111_2
% comb_C2_2 = [C2000_2_sq C2001_2_sq C2010_2_sq C2011_2_sq C2100_2_sq C2101_2_sq C2110_2_sq C2111_2_sq];


%000_2
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23000_2_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)

%001_2
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23001_2_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)


%010_2
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23010_2_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)

%011_2
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23011_2_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)

%100_2
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23100_2_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)

%101_2
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23101_2_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)


%110_2
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23110_2_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)


%111_2
fig = figure
set(fig, 'Position', [0 0 560 720])
Input = L23111_2_sq;
subplot(4,3,[1 2 4 5 7 8])
imagesc((1-Input))
colormap(gray)
subplot(4,3,[10 11])
bar(sum(Input(:,:),1))
xlim([1 tstop])
set(gca, 'XTick', [])
ylim([0 nDim*nDim*0.8])
subplot(4,3,[3 6 9])
firingCount = sum(Input(nDim*nDim:-1:1,:),2);
barh(firingCount)
xlim([1 13])
ylim([1 nDim*nDim])
set(gca, 'YTick', [])
tmpImg = zeros(nDim, nDim);
for y=1:nDim
    for x=1:nDim
        tmpImg(y,x) = firingCount((y-1)*nDim+x);
    end
end
subplot(4,3,12)
imagesc(tmpImg)







%combined
fig = figure
set(fig, 'Position', [0 0 560 720])
subplot(19,1,1:5)
Input = [comb_L5_2;comb_L23_2;comb_L4_2];
imagesc((1-Input(1:nDim*nDim,:)))
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylabel('V1_5 (900 cells)');

subplot(19,1,6)
bar(sum(Input(1:nDim*nDim,:),1))
xlim([0 (tstop+1)*8-1])
set(gca, 'XTick', [])

subplot(19,1,7:11)
imagesc((1-Input(nDim*nDim+1:nDim*nDim*2,:)))
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylabel('V1_2_/_3 (900 cells)');

subplot(19,1,12)
bar(sum(Input(nDim*nDim+1:nDim*nDim*2,:),1))
xlim([0 (tstop+1)*8-1])
set(gca, 'XTick', [])

subplot(19,1,13:17)
imagesc((1-Input(nDim*nDim*2+1:nDim*nDim*3,:)))
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylabel('V1_4 (900 cells)');

subplot(19,1,18)
bar(sum(Input(nDim*nDim*2+1:nDim*nDim*3,:),1))
xlim([0 (tstop+1)*8-1])
set(gca, 'XTick', [])

subplot(19,1,19)
image(img);
set(gca, 'XTick', [])
set(gca, 'YTick', [])
xlabel('input colors [240 ms/color] with 60ms intervals');











%combined 2
fig = figure

%000_on
tmp = [L5000_2_sq sum(L5000_2_sq(:,1:60),2)-sum(L5000_2_sq(:,140:190),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);

subplot(10,2,1:2:11)
imagesc((1-L5000_sort(1:30,:)));
colormap(gray)
set(gca, 'YTick', [])
set(gca, 'XTick', [])
ylim([1 30]);
subplot(10,2,13:2:17)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 20])
set(gca, 'XTick', [])
subplot(10,2,19)
image(img(1,1:5,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])

%000_off
tmp = [L5000_2_sq sum(L5000_2_sq(:,251:301),2).*1.3-sum(L5000_2_sq(:,110:160),2).*1.5 ];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);
subplot(10,2,2:2:12)
imagesc((1-L5000_sort(1:900,:)));
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylim([1 30]);
subplot(10,2,14:2:18)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 10])
set(gca, 'XTick', [])
subplot(10,2,20)
image(img(1,1:5,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])
%saveas(fig, 'L5000_on_off.png');
clf




%001_on
tmp = [L5001_2_sq sum(L5001_2_sq(:,1:60),2)-sum(L5001_2_sq(:,160:190),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);

subplot(10,2,1:2:11)
imagesc((1-L5000_sort(1:30,:)));
colormap(gray)
set(gca, 'YTick', [])
set(gca, 'XTick', [])
ylim([1 30]);
subplot(10,2,13:2:17)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 20])
set(gca, 'XTick', [])
subplot(10,2,19)
image(img(1,6:10,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])

%001_off
tmp = [L5001_2_sq sum(L5001_2_sq(:,251:301),2)-sum(L5001_2_sq(:,110:160),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);
subplot(10,2,2:2:12)
imagesc((1-L5000_sort(1:900,:)));
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylim([1 30]);
subplot(10,2,14:2:18)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 10])
set(gca, 'XTick', [])
subplot(10,2,20)
image(img(1,6:10,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])
clf


%010_on
tmp = [L5010_2_sq sum(L5010_2_sq(:,1:60),2)-sum(L5010_2_sq(:,250:300),2).*0.5];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);

subplot(10,2,1:2:11)
imagesc((1-L5000_sort(1:30,:)));
colormap(gray)
set(gca, 'YTick', [])
set(gca, 'XTick', [])
ylim([1 30]);
subplot(10,2,13:2:17)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 20])
set(gca, 'XTick', [])
subplot(10,2,19)
image(img(1,11:15,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])

%010_off
tmp = [L5010_2_sq sum(L5010_2_sq(:,251:301),2)-sum(L5010_2_sq(:,1:240),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);
subplot(10,2,2:2:12)
imagesc((1-L5000_sort(1:900,:)));
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylim([1 30]);
subplot(10,2,14:2:18)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 10])
set(gca, 'XTick', [])
subplot(10,2,20)
image(img(1,11:15,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])
clf


%011_on
tmp = [L5011_2_sq sum(L5011_2_sq(:,1:60),2)];%-sum(L5001_2_sq(:,160:190),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);

subplot(10,2,1:2:11)
imagesc((1-L5000_sort(1:30,:)));
colormap(gray)
set(gca, 'YTick', [])
set(gca, 'XTick', [])
ylim([1 30]);
subplot(10,2,13:2:17)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 20])
set(gca, 'XTick', [])
subplot(10,2,19)
image(img(1,16:20,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])

%011_off
tmp = [L5011_2_sq sum(L5011_2_sq(:,251:301),2)-sum(L5011_2_sq(:,1:240),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);
subplot(10,2,2:2:12)
imagesc((1-L5000_sort(1:900,:)));
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylim([1 30]);
subplot(10,2,14:2:18)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 10])
set(gca, 'XTick', [])
subplot(10,2,20)
image(img(1,16:20,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])
clf



%100_on
tmp = [L5100_2_sq sum(L5100_2_sq(:,1:60),2)-sum(L5100_2_sq(:,250:300),2).*0.5];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);

subplot(10,2,1:2:11)
imagesc((1-L5000_sort(1:30,:)));
colormap(gray)
set(gca, 'YTick', [])
set(gca, 'XTick', [])
ylim([1 30]);
subplot(10,2,13:2:17)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 20])
set(gca, 'XTick', [])
subplot(10,2,19)
image(img(1,21:25,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])

%100_off
tmp = [L5100_2_sq sum(L5100_2_sq(:,251:301),2)-sum(L5100_2_sq(:,1:240),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);
subplot(10,2,2:2:12)
imagesc((1-L5000_sort(1:900,:)));
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylim([1 30]);
subplot(10,2,14:2:18)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 10])
set(gca, 'XTick', [])
subplot(10,2,20)
image(img(1,21:25,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])
clf



%101_on
tmp = [L5101_2_sq sum(L5101_2_sq(:,1:60),2)-sum(L5101_2_sq(:,250:300),2).*0.5];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);

subplot(10,2,1:2:11)
imagesc((1-L5000_sort(1:30,:)));
colormap(gray)
set(gca, 'YTick', [])
set(gca, 'XTick', [])
ylim([1 30]);
subplot(10,2,13:2:17)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 20])
set(gca, 'XTick', [])
subplot(10,2,19)
image(img(1,26:30,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])

%101_off
tmp = [L5101_2_sq sum(L5101_2_sq(:,251:301),2)-sum(L5101_2_sq(:,1:240),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);
subplot(10,2,2:2:12)
imagesc((1-L5000_sort(1:900,:)));
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylim([1 30]);
subplot(10,2,14:2:18)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 10])
set(gca, 'XTick', [])
subplot(10,2,20)
image(img(1,26:30,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])
clf


%110_on
tmp = [L5110_2_sq sum(L5110_2_sq(:,1:60),2)-sum(L5110_2_sq(:,250:300),2).*0.5];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);

subplot(10,2,1:2:11)
imagesc((1-L5000_sort(1:30,:)));
colormap(gray)
set(gca, 'YTick', [])
set(gca, 'XTick', [])
ylim([1 30]);
subplot(10,2,13:2:17)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 20])
set(gca, 'XTick', [])
subplot(10,2,19)
image(img(1,31:35,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])

%101_off
tmp = [L5110_2_sq sum(L5110_2_sq(:,251:301),2)-sum(L5110_2_sq(:,1:240),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);
subplot(10,2,2:2:12)
imagesc((1-L5000_sort(1:900,:)));
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylim([1 30]);
subplot(10,2,14:2:18)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 10])
set(gca, 'XTick', [])
subplot(10,2,20)
image(img(1,31:35,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])
clf


%111_on
tmp = [L5111_2_sq sum(L5111_2_sq(:,1:60),2)-sum(L5111_2_sq(:,250:300),2).*0.5];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);

subplot(10,2,1:2:11)
imagesc((1-L5000_sort(1:30,:)));
colormap(gray)
set(gca, 'YTick', [])
set(gca, 'XTick', [])
ylim([1 30]);
subplot(10,2,13:2:17)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 20])
set(gca, 'XTick', [])
subplot(10,2,19)
image(img(1,36:40,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])

%101_off
tmp = [L5111_2_sq sum(L5111_2_sq(:,251:301),2)-sum(L5111_2_sq(:,1:240),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);
subplot(10,2,2:2:12)
imagesc((1-L5000_sort(1:900,:)));
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylim([1 30]);
subplot(10,2,14:2:18)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 10])
set(gca, 'XTick', [])
subplot(10,2,20)
image(img(1,36:40,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])
clf



fig = figure
%001 and 110 comb
tmp = [L5110_2_sq L5001_2_sq sum(L5110_2_sq(:,1:60),2)-sum(L5001_2_sq(:,110:160),2)-sum(L5110_2_sq(:,250:300),2)];
L5110_sort = sortrows(tmp,603);
L5110_sort = L5110_sort(601:900,1:602);
tmp = [L5110_sort sum(L5110_sort(:,551:602),2)-sum(L5110_sort(:,301:540),2)];
L5110_sort = sortrows(tmp,603);
L5110_sort = L5110_sort(271:300,1:602);
subplot(10,1,1:6)
imagesc((1-L5110_sort(:,:)));
colormap(gray)
set(gca, 'YTick', [])
set(gca, 'XTick', [])
subplot(10,1,7:9)
bar(sum(L5110_sort(:,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 601]);
set(gca, 'XTick', [])
subplot(10,1,10)
image([img(1,31:35,:) img(1,6:10,:)]);
set(gca, 'XTick', [])
set(gca, 'YTick', [])
saveas(fig, 'L5110_on_L5001_off.png');
clf


%100 and 011 comb
tmp = [L5100_2_sq L5010_2_sq sum(L5100_2_sq(:,1:60),2)-sum(L5100_2_sq(:,250:300),2)];
L5100_sort = sortrows(tmp,603);
L5100_sort = L5100_sort(601:900,1:602);
tmp = [L5100_sort sum(L5100_sort(:,551:602),2)-sum(L5100_sort(:,301:540),2)];
L5100_sort = sortrows(tmp,603);
L5100_sort = L5100_sort(271:300,1:602);

subplot(10,1,1:6)
imagesc((1-L5100_sort(:,:)));
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
subplot(10,1,7:9)
bar(sum(L5100_sort(:,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 601]);
%ylim([0 20])
set(gca, 'XTick', [])
subplot(10,1,10)
image([img(1,21:25,:) img(1,11:15,:)]);
set(gca, 'XTick', [])
set(gca, 'YTick', [])
saveas(fig, 'L5100_on_L5011_off.png');
clf






%001_off
tmp = [L5001_2_sq sum(L5001_2_sq(:,251:301),2)-sum(L5001_2_sq(:,110:160),2)];
L5000_sort = -sortrows(-tmp,302);
L5000_sort = L5000_sort(:,1:301);
subplot(10,2,2:2:12)
imagesc((1-L5000_sort(1:900,:)));
colormap(gray)
set(gca, 'XTick', [])
set(gca, 'YTick', [])
ylim([1 30]);
subplot(10,2,14:2:18)
bar(sum(L5000_sort(1:30,:),1))
%plot(sum(L5000_sort(800:900,:),1))
xlim([0 300]);
%ylim([0 10])
set(gca, 'XTick', [])
subplot(10,2,20)
image(img(1,6:10,:));
set(gca, 'XTick', [])
set(gca, 'YTick', [])
clf







% fig = figure
% set(fig, 'Position', [0 0 560 720])
% subplot(19,1,1:6)
% Input = comb_L_2;
% imagesc((1-Input))
% colormap(gray)
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])
% ylabel('L channel (R+G)')
% 
% subplot(19,1,7:12)
% Input = comb_C1_2;
% imagesc((1-Input))
% colormap(gray)
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])
% ylabel('C1 channel (R-G)')
% 
% subplot(19,1,13:18)
% Input = comb_C2_2;
% imagesc((1-Input))
% colormap(gray)
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])
% ylabel('C2 channel (R+G-B)')
% 
% 
% subplot(19,1,19)
% image(img);
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])
% xlabel('input colors [240 ms/color] with 60ms intervals')


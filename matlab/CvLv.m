nDim = 30
tstop = 1000
dt = 0.025

L23000 = load('L23_firing_1000_[0, 0, 0].txt');
L23001 = load('L23_firing_1000_[0, 0, 1].txt');
L23010 = load('L23_firing_1000_[0, 1, 0].txt');
L23011 = load('L23_firing_1000_[0, 1, 1].txt');
L23100 = load('L23_firing_1000_[1, 0, 0].txt');
L23101 = load('L23_firing_1000_[1, 0, 1].txt');
L23110 = load('L23_firing_1000_[1, 1, 0].txt');
L23111 = load('L23_firing_1000_[1, 1, 1].txt');
L23000_sq = zeros(nDim*nDim,tstop+1);
L23001_sq = zeros(nDim*nDim,tstop+1);
L23010_sq = zeros(nDim*nDim,tstop+1);
L23011_sq = zeros(nDim*nDim,tstop+1);
L23100_sq = zeros(nDim*nDim,tstop+1);
L23101_sq = zeros(nDim*nDim,tstop+1);
L23110_sq = zeros(nDim*nDim,tstop+1);
L23111_sq = zeros(nDim*nDim,tstop+1);
for t = 0:tstop-1
    index = 2+(t)*40;
    L23000_sq(:,t+1) = mean(L23000(:,index:index+40),2);
    L23001_sq(:,t+1) = mean(L23001(:,index:index+40),2);
    L23010_sq(:,t+1) = mean(L23010(:,index:index+40),2);
    L23011_sq(:,t+1) = mean(L23011(:,index:index+40),2);
    L23100_sq(:,t+1) = mean(L23100(:,index:index+40),2);
    L23101_sq(:,t+1) = mean(L23101(:,index:index+40),2);
    L23110_sq(:,t+1) = mean(L23110(:,index:index+40),2);
    L23111_sq(:,t+1) = mean(L23111(:,index:index+40),2);
end

L23000 = zeros(1);
L23001 = zeros(1);
L23010 = zeros(1);
L23011 = zeros(1);
L23100 = zeros(1);
L23101 = zeros(1);
L23110 = zeros(1);
L23111 = zeros(1);

nSegs = 10;

%000
data = L23100_sq;
ISI_100 = zeros(nDim*nDim,nSegs);
for index=1:nDim*nDim
    prevt = 0;
    iCount = 1;
    for t = 1:tstop+1
        if ISI_100(index,nSegs)~=0
            break;
        elseif data(index,t)>0
            if (prevt==0)
                prevt=t;
            elseif(t-prevt>5)
                ISI_100(index,iCount) = (t-prevt);
                iCount = iCount+1;
                prevt=t;
            end
        end
    end
end


Tmean = zeros(nDim*nDim,1);
Tstd = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    Tmean(index) = mean(ISI_100(index,find(ISI_100(index,:)~=0)));
    Tstd(index) = std(ISI_100(index,find(ISI_100(index,:)~=0)));
end
%Tstd(find(isnan(Tstd))) = 1;
%Tmean(find(isnan(Tmean))) = -1;
Cv = Tstd./Tmean
%tmp = sortrows([Cv data]);
%tmp2 = tmp(find(tmp(:,1)>1),2:tstop+2);
tmp = data(find(tmp(:,1)>1),:);
tmp2 = data(find(tmp(:,1)<=1),:);

figure;
subplot(8,2,1:2:9)
imagesc(tmp(:,:))
subplot(8,2,11:2:15)
bar(sum(tmp(:,:),1))
xlim([0 (tstop+1)])
%ylim([-0.1 10])

subplot(8,2,2:2:10)
imagesc(tmp2(:,:))
subplot(8,2,12:2:16)
bar(sum(tmp2(:,:),1))
xlim([0 (tstop+1)])



Lv = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    tot = 0;
    for i=1:nSegs-1
        if(power(ISI_100(index,i)+ISI_100(index,i+1),2)~=0)
            tot = tot + 3*power(ISI_100(index,i)-ISI_100(index,i+1),2)/power(ISI_100(index,i)+ISI_100(index,i+1),2);
        end
    end
    Lv(index) = 1/(nSegs-1)*tot;
end

tmp = sortrows([Lv data]);
figure;
imagesc(tmp(:,2:tstop+2))




CvAv = mean(Cv(find(~isnan(Cv))));
LvAvg = mean(Lv(find(~isnan(Lv))));
[CvAv LvAvg]

%001
data = L23001_sq;
ISI_100 = zeros(nDim*nDim,nSegs);
for index=1:nDim*nDim
    prevt = 0;
    iCount = 1;
    for t = 1:tstop+1
        if ISI_100(index,nSegs)~=0
            break;
        elseif data(index,t)>0
            if (prevt==0)
                prevt=t;
            elseif(t-prevt>5)
                ISI_100(index,iCount) = (t-prevt);
                iCount = iCount+1;
                prevt=t;
            end
        end
    end
end


Tmean = zeros(nDim*nDim,1);
Tstd = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    Tmean(index) = mean(ISI_100(index,find(ISI_100(index,:)~=0)));
    Tstd(index) = std(ISI_100(index,find(ISI_100(index,:)~=0)));
end
Tstd(find(isnan(Tstd))) = 1;
Tmean(find(isnan(Tmean))) = -1;
Cv = Tstd./Tmean;


Lv = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    tot = 0;
    for i=1:nSegs-1
        tot = tot + 3*power(ISI_100(index,i)-ISI_100(index,i+1),2)/power(ISI_100(index,i)+ISI_100(index,i+1),2);
    end
    Lv(index) = 1/(nSegs-1)*tot;
end
CvAv = mean(Cv(find(~isnan(Cv))));
LvAvg = mean(Lv(find(~isnan(Lv))));
[CvAv LvAvg]


%010
data = L23010_sq;
ISI_100 = zeros(nDim*nDim,nSegs);
for index=1:nDim*nDim
    prevt = 0;
    iCount = 1;
    for t = 1:tstop+1
        if ISI_100(index,nSegs)~=0
            break;
        elseif data(index,t)>0
            if (prevt==0)
                prevt=t;
            elseif(t-prevt>5)
                ISI_100(index,iCount) = (t-prevt);
                iCount = iCount+1;
                prevt=t;
            end
        end
    end
end


Tmean = zeros(nDim*nDim,1);
Tstd = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    Tmean(index) = mean(ISI_100(index,find(ISI_100(index,:)~=0)));
    Tstd(index) = std(ISI_100(index,find(ISI_100(index,:)~=0)));
end
Tstd(find(isnan(Tstd))) = 1;
Tmean(find(isnan(Tmean))) = -1;
Cv = Tstd./Tmean;


Lv = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    tot = 0;
    for i=1:nSegs-1
        tot = tot + 3*power(ISI_100(index,i)-ISI_100(index,i+1),2)/power(ISI_100(index,i)+ISI_100(index,i+1),2);
    end
    Lv(index) = 1/(nSegs-1)*tot;
end
CvAv = mean(Cv(find(~isnan(Cv))));
LvAvg = mean(Lv(find(~isnan(Lv))));
[CvAv LvAvg]


%011
data = L23011_sq;
ISI_100 = zeros(nDim*nDim,nSegs);
for index=1:nDim*nDim
    prevt = 0;
    iCount = 1;
    for t = 1:tstop+1
        if ISI_100(index,nSegs)~=0
            break;
        elseif data(index,t)>0
            if (prevt==0)
                prevt=t;
            elseif(t-prevt>5)
                ISI_100(index,iCount) = (t-prevt);
                iCount = iCount+1;
                prevt=t;
            end
        end
    end
end


Tmean = zeros(nDim*nDim,1);
Tstd = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    Tmean(index) = mean(ISI_100(index,find(ISI_100(index,:)~=0)));
    Tstd(index) = std(ISI_100(index,find(ISI_100(index,:)~=0)));
end
Tstd(find(isnan(Tstd))) = 1;
Tmean(find(isnan(Tmean))) = -1;
Cv = Tstd./Tmean;


Lv = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    tot = 0;
    for i=1:nSegs-1
        tot = tot + 3*power(ISI_100(index,i)-ISI_100(index,i+1),2)/power(ISI_100(index,i)+ISI_100(index,i+1),2);
    end
    Lv(index) = 1/(nSegs-1)*tot;
end
CvAv = mean(Cv(find(~isnan(Cv))));
LvAvg = mean(Lv(find(~isnan(Lv))));
[CvAv LvAvg]


%100
data = L23100_sq;
ISI_100 = zeros(nDim*nDim,nSegs);
for index=1:nDim*nDim
    prevt = 0;
    iCount = 1;
    for t = 1:tstop+1
        if ISI_100(index,nSegs)~=0
            break;
        elseif data(index,t)>0
            if (prevt==0)
                prevt=t;
            elseif(t-prevt>5)
                ISI_100(index,iCount) = (t-prevt);
                iCount = iCount+1;
                prevt=t;
            end
        end
    end
end


Tmean = zeros(nDim*nDim,1);
Tstd = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    Tmean(index) = mean(ISI_100(index,find(ISI_100(index,:)~=0)));
    Tstd(index) = std(ISI_100(index,find(ISI_100(index,:)~=0)));
end
Tstd(find(isnan(Tstd))) = 1;
Tmean(find(isnan(Tmean))) = -1;
Cv = Tstd./Tmean;


Lv = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    tot = 0;
    for i=1:nSegs-1
        tot = tot + 3*power(ISI_100(index,i)-ISI_100(index,i+1),2)/power(ISI_100(index,i)+ISI_100(index,i+1),2);
    end
    Lv(index) = 1/(nSegs-1)*tot;
end
CvAv = mean(Cv(find(~isnan(Cv))));
LvAvg = mean(Lv(find(~isnan(Lv))));
[CvAv LvAvg]



%101
data = L23101_sq;
ISI_100 = zeros(nDim*nDim,nSegs);
for index=1:nDim*nDim
    prevt = 0;
    iCount = 1;
    for t = 1:tstop+1
        if ISI_100(index,nSegs)~=0
            break;
        elseif data(index,t)>0
            if (prevt==0)
                prevt=t;
            elseif(t-prevt>5)
                ISI_100(index,iCount) = (t-prevt);
                iCount = iCount+1;
                prevt=t;
            end
        end
    end
end


Tmean = zeros(nDim*nDim,1);
Tstd = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    Tmean(index) = mean(ISI_100(index,find(ISI_100(index,:)~=0)));
    Tstd(index) = std(ISI_100(index,find(ISI_100(index,:)~=0)));
end
Tstd(find(isnan(Tstd))) = 1;
Tmean(find(isnan(Tmean))) = -1;
Cv = Tstd./Tmean;


Lv = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    tot = 0;
    for i=1:nSegs-1
        tot = tot + 3*power(ISI_100(index,i)-ISI_100(index,i+1),2)/power(ISI_100(index,i)+ISI_100(index,i+1),2);
    end
    Lv(index) = 1/(nSegs-1)*tot;
end
CvAv = mean(Cv(find(~isnan(Cv))));
LvAvg = mean(Lv(find(~isnan(Lv))));
[CvAv LvAvg]



%110
data = L23110_sq;
ISI_100 = zeros(nDim*nDim,nSegs);
for index=1:nDim*nDim
    prevt = 0;
    iCount = 1;
    for t = 1:tstop+1
        if ISI_100(index,nSegs)~=0
            break;
        elseif data(index,t)>0
            if (prevt==0)
                prevt=t;
            elseif(t-prevt>5)
                ISI_100(index,iCount) = (t-prevt);
                iCount = iCount+1;
                prevt=t;
            end
        end
    end
end


Tmean = zeros(nDim*nDim,1);
Tstd = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    Tmean(index) = mean(ISI_100(index,find(ISI_100(index,:)~=0)));
    Tstd(index) = std(ISI_100(index,find(ISI_100(index,:)~=0)));
end
Tstd(find(isnan(Tstd))) = 1;
Tmean(find(isnan(Tmean))) = -1;
Cv = Tstd./Tmean;


Lv = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    tot = 0;
    for i=1:nSegs-1
        tot = tot + 3*power(ISI_100(index,i)-ISI_100(index,i+1),2)/power(ISI_100(index,i)+ISI_100(index,i+1),2);
    end
    Lv(index) = 1/(nSegs-1)*tot;
end
CvAv = mean(Cv(find(~isnan(Cv))));
LvAvg = mean(Lv(find(~isnan(Lv))));
[CvAv LvAvg]


%111
data = L23111_sq;
ISI_100 = zeros(nDim*nDim,nSegs);
for index=1:nDim*nDim
    prevt = 0;
    iCount = 1;
    for t = 1:tstop+1
        if ISI_100(index,nSegs)~=0
            break;
        elseif data(index,t)>0
            if (prevt==0)
                prevt=t;
            elseif(t-prevt>5)
                ISI_100(index,iCount) = (t-prevt);
                iCount = iCount+1;
                prevt=t;
            end
        end
    end
end


Tmean = zeros(nDim*nDim,1);
Tstd = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    Tmean(index) = mean(ISI_100(index,find(ISI_100(index,:)~=0)));
    Tstd(index) = std(ISI_100(index,find(ISI_100(index,:)~=0)));
end
Tstd(find(isnan(Tstd))) = 1;
Tmean(find(isnan(Tmean))) = -1;
Cv = Tstd./Tmean;


Lv = zeros(nDim*nDim,1);
for index=1:nDim*nDim
    tot = 0;
    for i=1:nSegs-1
        tot = tot + 3*power(ISI_100(index,i)-ISI_100(index,i+1),2)/power(ISI_100(index,i)+ISI_100(index,i+1),2);
    end
    Lv(index) = 1/(nSegs-1)*tot;
end
CvAv = mean(Cv(find(~isnan(Cv))));
LvAvg = mean(Lv(find(~isnan(Lv))));
[CvAv LvAvg]
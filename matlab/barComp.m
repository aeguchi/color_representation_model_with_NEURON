FR1 = load('L23_FR_1000_[0, 0, 0].txt');
FR2 = load('L23_FR_1000_[0, 0, 1].txt');
FR3 = load('L23_FR_1000_[0, 1, 0].txt');
FR4 = load('L23_FR_1000_[0, 1, 1].txt');
FR5 = load('L23_FR_1000_[1, 0, 0].txt');
FR6 = load('L23_FR_1000_[1, 0, 1].txt');
FR7 = load('L23_FR_1000_[1, 1, 0].txt');
FR8 = load('L23_FR_1000_[1, 1, 1].txt');

% FR1 = load('L23_FR_0_[0, 0, 0].txt');
% FR2 = load('L23_FR_0_[0, 0, 1].txt');
% FR3 = load('L23_FR_0_[0, 1, 0].txt');
% FR4 = load('L23_FR_0_[0, 1, 1].txt');
% FR5 = load('L23_FR_0_[1, 0, 0].txt');
% FR6 = load('L23_FR_0_[1, 0, 1].txt');
% FR7 = load('L23_FR_0_[1, 1, 0].txt');
% FR8 = load('L23_FR_0_[1, 1, 1].txt');


FRavg = (FR1+FR2+FR3+FR4+FR5+FR6+FR7+FR8)./8;


FR = [FR1; FR2; FR3; FR4; FR5; FR6; FR7; FR8];
%num_transforms, num_objects, cell_x, cell_y

num_cells = sqrt(length(FR))
%num_bins = max(FR(:))+1
max_FR = max(FR(:))
num_bins = 10
num_stimulus = length(FR(:,1))/length(FR1(:,1))
num_transforms = length(FR1(:,1))

if (max_FR==0)
    return
end

binMatrix = zeros(num_cells,num_cells,num_stimulus,num_bins);
sumPerObj = num_transforms
sumPerCell = num_stimulus*num_transforms


for y = 1:num_cells
    for x=1:num_cells
        index = y*(num_cells-1)+x;
        for stim = 1:num_stimulus
            for trans = 1:num_transforms
                binMatrix(x,y,stim,int8(FR((stim-1)*num_transforms+trans,index)/max_FR*(num_bins-1))+1)=binMatrix(x,y,stim,int8(FR((stim-1)*num_transforms+trans,index)/max_FR*(num_bins-1))+1)+1;
                %test = [int8(FR(color,index)/max_FR*(num_bins+1))+1]
            end
        end
    end
end


sumPerBin = zeros(num_cells,num_cells,num_bins);
IRs = zeros(num_cells,num_cells,num_stimulus);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single-cell information analysis      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through all cells to calculate single cell information
for cell_x=1:num_cells
    for cell_y=1:num_cells
        % For each cell, count the number of transforms per bin
        for bin=1:num_bins
            sumPerBin(cell_x,cell_y,bin)=sum(binMatrix(cell_x,cell_y,:,bin));
        end


        % Calculate the information for cell_x cell_y per stimulus
        for stimulus=1:num_stimulus
            for bin=1:num_bins
                Pr = sumPerBin(cell_x,cell_y,bin)/sumPerCell;
                Prs = binMatrix(cell_x,cell_y,stimulus,bin)/sumPerObj;
                if(Pr~=0&&Prs~=0)
                    IRs(cell_x,cell_y,stimulus)=IRs(cell_x,cell_y,stimulus)+(Prs*(log2(Prs/Pr)));
                    %[Pr Prs Prs*(log2(Prs/Pr))]
                end
            end
        end
    end
end
h=figure;
powVar = 1;
tmpImg = zeros(num_cells, num_cells,3);
for y=1:num_cells
    for x = 1:num_cells
        plotIndex = mod((y-1)*num_cells+x-1,64)+1;
        subplot(8,8,plotIndex)
%         xlim([1 8]);
%         ylim([0 3]);
        maxVal = max(IRs(y,x,:));
        if(maxVal>0)
            %clf
            for color = 1:8
                tmpVar = zeros(1,8);
                if(color==1)
                    tmpVar(color) = IRs(x,y,color);
                    h1 = bar(tmpVar);
                    set(h1,'facecolor',[0 0 0]) 
%                     tmpVar(color) = FRavg(color,(y-1)*num_cells+x)-mean(FRavg(:,(y-1)*num_cells+x));
%                     stem(tmpVar)
                elseif(color==2)
                    tmpVar(color) = IRs(x,y,color);
                    h1 = bar(tmpVar);
                    set(h1,'facecolor',[0 0 1]) 
%                     tmpVar(color) = FRavg(color,(y-1)*num_cells+x)-mean(FRavg(:,(y-1)*num_cells+x));
%                     stem(tmpVar)
                elseif(color==3)
                    tmpVar(color) = IRs(x,y,color);
                    h1 = bar(tmpVar);
                    set(h1,'facecolor',[0 1 0]) 
%                     tmpVar(color) = FRavg(color,(y-1)*num_cells+x)-mean(FRavg(:,(y-1)*num_cells+x));
%                     stem(tmpVar)
                elseif(color==4)
                    tmpVar(color) = IRs(x,y,color);
                    h1 = bar(tmpVar);
                    set(h1,'facecolor',[0 1 1]) 
%                     tmpVar(color) = FRavg(color,(y-1)*num_cells+x)-mean(FRavg(:,(y-1)*num_cells+x));
%                     stem(tmpVar)
                elseif(color==5)
                    tmpVar(color) = IRs(x,y,color);
                    h1 = bar(tmpVar);
                    set(h1,'facecolor',[1 0 0])  
%                     tmpVar(color) = FRavg(color,(y-1)*num_cells+x)-mean(FRavg(:,(y-1)*num_cells+x));
%                     stem(tmpVar)               
                elseif(color==6)
                    tmpVar(color) = IRs(x,y,color);
                    h1 = bar(tmpVar);
                    set(h1,'facecolor',[1 0 1]) 
%                     tmpVar(color) = FRavg(color,(y-1)*num_cells+x)-mean(FRavg(:,(y-1)*num_cells+x));
%                     stem(tmpVar)
                elseif(color==7)
                    tmpVar(color) = IRs(x,y,color);
                    h1 = bar(tmpVar);
                    set(h1,'facecolor',[1 1 0]) 
%                     tmpVar(color) = FRavg(color,(y-1)*num_cells+x)-mean(FRavg(:,(y-1)*num_cells+x));
%                     stem(tmpVar)
                elseif(color==8)
                    tmpVar(color) = IRs(x,y,color);
                    h1 = bar(tmpVar);
                    set(h1,'facecolor',[1 1 1]) 
%                     tmpVar(color) = FRavg(color,(y-1)*num_cells+x)-mean(FRavg(:,(y-1)*num_cells+x));
%                     stem(tmpVar)
                end
                hold on
            end
            xlim([0 9])
            ylim([0 3])
        else
            tmpImg(y,x,:) = [0.5 0.5 0.5];
        end
        if(plotIndex==64)
            saveas(h,['infoBar_' num2str(y) '_' num2str(x) '.png']) 
            clf;
%             figure
        end
    end
end
figure
image(tmpImg)
        


% fig = figure
% %set(fig, 'Position', [0 0 560 720])
% Input = L23000_2_sq;
% firingCount = sum(Input(nDim*nDim:-1:1,:),2);
% Input = L23001_2_sq;
% firingCount = firingCount + sum(Input(nDim*nDim:-1:1,:),2);
% Input = L23010_2_sq;
% firingCount = firingCount + sum(Input(nDim*nDim:-1:1,:),2);
% Input = L23011_2_sq;
% firingCount = firingCount + sum(Input(nDim*nDim:-1:1,:),2);
% Input = L23100_2_sq;
% firingCount = firingCount + sum(Input(nDim*nDim:-1:1,:),2);
% Input = L23101_2_sq;
% firingCount = firingCount + sum(Input(nDim*nDim:-1:1,:),2);
% Input = L23110_2_sq;
% firingCount = firingCount + sum(Input(nDim*nDim:-1:1,:),2);
% Input = L23111_2_sq;
% firingCount = firingCount + sum(Input(nDim*nDim:-1:1,:),2);
% 
% tmpImg = zeros(nDim, nDim,3);
% for y=1:nDim
%     for x=1:nDim
%         tmpImg(y,x,:) = [firingCount((y-1)*nDim+x) firingCount((y-1)*nDim+x) firingCount((y-1)*nDim+x)];
%     end
% end
% tmpImg = tmpImg/max(tmpImg(:));
% imagesc(tmpImg)
% figure
% 
% for i = 1:15
%     index = topIRs(i,1);
%     y = ceil(index/30);
%     x = mod(index-1,30)+1;
%     tmpImg(y,x,:)=[0 0 0];
% end
% 
% for i = 1:15
%     index = topIRs(i,2);
%     y = ceil(index/30);
%     x = mod(index-1,30)+1;
%     tmpImg(y,x,:)=[0 0 1];
% end
% 
% for i = 1:15
%     index = topIRs(i,3);
%     y = ceil(index/30);
%     x = mod(index-1,30)+1;
%     tmpImg(y,x,:)=[0 1 0];
% end
% 
% for i = 1:15
%     index = topIRs(i,4);
%     y = ceil(index/30);
%     x = mod(index-1,30)+1;
%     tmpImg(y,x,:)=[0 1 1];
% end
% 
% for i = 1:15
%     index = topIRs(i,5);
%     y = ceil(index/30);
%     x = mod(index-1,30)+1;
%     tmpImg(y,x,:)=[1 0 0];
% end
% 
% for i = 1:15
%     index = topIRs(i,6);
%     y = ceil(index/30);
%     x = mod(index-1,30)+1;
%     tmpImg(y,x,:)=[1 0 1];
% end
% 
% for i = 1:15
%     index = topIRs(i,7);
%     y = ceil(index/30);
%     x = mod(index-1,30)+1;
%     tmpImg(y,x,:)=[1 1 0];
% end
% 
% for i = 1:15
%     index = topIRs(i,1);
%     y = ceil(index/30);
%     x = mod(index-1,30)+1;
%     tmpImg(y,x,:)=[1 1 1];
% end
% 
% imagesc(tmpImg)

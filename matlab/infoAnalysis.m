function infoAnalysis(fileNames,nc_max)
%USAGE:
% infoAnalysis({'0' '1000'}, nc_max);

    multi_cell_analysis = 1
    figure;
    
    tmp = load(['L5_FR_',num2str(fileNames{1}),'_[0, 0, 0].txt']);
    num_cells = sqrt(length(tmp));
    nData = length(fileNames);
    SingleCellPlot = zeros(nData,num_cells,num_cells);
    MultiCellPlotP = zeros(nData,nc_max+1);
    MultiCellPlotQ = zeros(nData,nc_max+1);
    
    for dataIndex = 1: length(fileNames)

        if(dataIndex~=3)
            %build binMatrix(cell_x,cell_y,:,bin) x y trans bin
            FR1 = load(['L5_FR_',num2str(fileNames{dataIndex}),'_[0, 0, 0].txt']);
            FR2 = load(['L5_FR_',num2str(fileNames{dataIndex}),'_[0, 0, 1].txt']);
            FR3 = load(['L5_FR_',num2str(fileNames{dataIndex}),'_[0, 1, 0].txt']);
            FR4 = load(['L5_FR_',num2str(fileNames{dataIndex}),'_[0, 1, 1].txt']);
            FR5 = load(['L5_FR_',num2str(fileNames{dataIndex}),'_[1, 0, 0].txt']);
            FR6 = load(['L5_FR_',num2str(fileNames{dataIndex}),'_[1, 0, 1].txt']);
            FR7 = load(['L5_FR_',num2str(fileNames{dataIndex}),'_[1, 1, 0].txt']);
            FR8 = load(['L5_FR_',num2str(fileNames{dataIndex}),'_[1, 1, 1].txt']);
        else
            FR1 = load(['L4_FR_',num2str(fileNames{dataIndex}),'_[0, 0, 0].txt']);
            FR2 = load(['L4_FR_',num2str(fileNames{dataIndex}),'_[0, 0, 1].txt']);
            FR3 = load(['L4_FR_',num2str(fileNames{dataIndex}),'_[0, 1, 0].txt']);
            FR4 = load(['L4_FR_',num2str(fileNames{dataIndex}),'_[0, 1, 1].txt']);
            FR5 = load(['L4_FR_',num2str(fileNames{dataIndex}),'_[1, 0, 0].txt']);
            FR6 = load(['L4_FR_',num2str(fileNames{dataIndex}),'_[1, 0, 1].txt']);
            FR7 = load(['L4_FR_',num2str(fileNames{dataIndex}),'_[1, 1, 0].txt']);
            FR8 = load(['L4_FR_',num2str(fileNames{dataIndex}),'_[1, 1, 1].txt']);
        end

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

        % Order by information content, descending
        IRs_tmp = sort(reshape(max(IRs,[],3),1,num_cells*num_cells), 'descend');%find max IRs for each


        SingleCellPlot(dataIndex,:) = IRs_tmp;


        firingRates = zeros(num_transforms, num_stimulus, num_cells, num_cells);
        for trans = 1:num_transforms
            for stim = 1:num_stimulus
                for y = 1:num_cells
                    for x=1:num_cells
                        firingRates(trans,stim,x,y) = FR((stim-1)*num_transforms+trans,(y-1)*num_cells+x)/max_FR;
                    end
                end
            end
        end

        randomPickForMulti=0
        UseMeanForMulti = 1
        pq_r = zeros(num_stimulus);  %prob(s') temporally used in the decoing process
        Ps = 1/num_stimulus; %Prob(s)
        Iss2_Q_avg = zeros(nc_max); %average quantized info for n cells; I(s,s')
        Iss2_P_avg = zeros(nc_max);% average smoothed info for n cells; I(s,s')


        num_samples = 5;    %num samples from Max-Cells of each stimulus for multi-cell info. analysis
        sampleMulti = 1;    %default is 1; this is to increase the sampling pool size to deal with larger layers.
        %nc_max = 20;%num_samples*num_stimulus;   % max ensemble size for multi-cell info. analysis


        IRs_topC = zeros(num_samples*sampleMulti,num_stimulus);  %to store top (num_samples) of cell_no's for each object.

        if(multi_cell_analysis)
            disp(['** multiple-cell information analysis **']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % multiple-cell information analysis    %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



            IRs_reshaped = [reshape([1:1:(num_cells*num_cells)],num_cells*num_cells,1) reshape(IRs,num_cells*num_cells,num_stimulus)]; %add index to represent cell_no
            for obj=1:num_stimulus %sort the IRs table and build a table of cell_no according to the amount of single cell info.
                IRs_sorted = sortrows(IRs_reshaped,-(obj+1));

                if(randomPickForMulti==0)
                    %            cellsWithMaxInf = find(IRs_sorted==IRs_sorted(1));
                    %            if length(cellsWithMaxInf)>=num_samples*sampleMulti
                    %                IRs_topC(:,obj) = IRs_sorted(cellsWithMaxInf(randperm(num_cells*num_cells,num_samples*sampleMulti)));
                    %            else
                    IRs_topC(:,obj) = IRs_sorted(1:num_samples*sampleMulti);
                    %            end
                elseif(randomPickForMulti==1)%just for testing purpose
                    IRs_topC(:,obj) = IRs_reshaped(randperm(num_cells*num_cells,num_samples*sampleMulti));
                end

            end

            IRs_topC_lined = reshape(IRs_topC(:,:),num_samples*sampleMulti*num_stimulus,1);

            if(UseMeanForMulti == 0)
                dp = 1/(num_stimulus * num_transforms);
            else
                dp = 1/num_stimulus;
            end
            for nc=1:nc_max

                disp([num2str(nc) '/' num2str(nc_max)]);

                %niter = 100*(nc_max/2);
                niter = 100*(nc_max - nc + 1);

                Iss2_Q_avg(nc) = 0.;
                Iss2_P_avg(nc) = 0.;


                testAvgP = zeros(num_stimulus,num_stimulus);
                testAvgQ = zeros(num_stimulus,num_stimulus);

                for iter = 1:niter
                    e = IRs_topC_lined(randperm(num_samples*sampleMulti*num_stimulus,nc)); %randomly pick nc number of cells which have max IRs
                    e2(:).cell_x = mod(e(:)-1,num_cells)+1;
                    e2(:).cell_y = floor((e(:)-1)/num_cells)+1;

                    %reset probability tables and info values
                    Iss2_Q = 0.;     %quantized raw info (< frequencies)
                    Iss2_P = 0.;  %smoothed raw info (< probabilities)

                    Pq = zeros(num_stimulus);%mean assigned probability P(s') (smoothed)
                    Psq = zeros(num_stimulus,num_stimulus);%probability table P(s,s')

                    Qq  = zeros(num_stimulus);%frequency of each predicted s Q(s') (extract the best)
                    Qsq = zeros(num_stimulus,num_stimulus);%frequency table Q(s,s')

                    ra = zeros(nc, num_stimulus);%(training set) averages
                    sd = zeros(nc, num_stimulus);%(training set) variances

                    if (UseMeanForMulti == 0)
                        rc = zeros(nc, num_stimulus);%current response
                        for s=1:num_stimulus
                            for tt=1:num_transforms
                                for ss=1:num_stimulus
                                    rc(:,ss) = transpose(reshape(firingRates(tt,ss,e(:)),1,nc));%get currecnt firing rates for each cell when exposed to object s at transform of tt
                                    ra(:,ss) = transpose(reshape(mean(firingRates(find([1:num_transforms]~=tt),ss,e(:))),1,nc));%get average firing rates of each cell when exposed to object s at other transforms
                                    sd(:,ss) = transpose(reshape(std(firingRates(find([1:num_transforms]~=tt),ss,e(:))),1,nc));%get standard deviations of the fr at other transforms
                                end

                                decode();

                                %dp is P(rc(s,t)) = 1 / (num_stimulus * num_transforms)
                                %P(s,s')= P(s'|rc(s,t)) * P(rc(s,t))

                                Psq(s,:) = Psq(s,:)+dp*reshape(pq_r(1:num_stimulus),1,num_stimulus);%probability table
                                Pq(:) = Pq(:) + dp * pq_r(:);%mean assigned probability
                                Qsq(s,best_s) = Qsq(s,best_s)+dp;%frequency table P(s,s')
                                Qq(best_s) = Qq(best_s) + dp;%frequency of each predicted s P(s')

                            end
                        end
                    else
                        for ss=1:num_stimulus
                            ra(:,ss) = transpose(reshape(mean(firingRates(:,ss,e(:))),1,nc));%get average firing rates of each cell when exposed to object s at other transforms
                            sd(:,ss) = transpose(reshape(std(firingRates(:,ss,e(:))),1,nc));%get standard deviations of the fr at other transforms
                        end



                        for s=1:num_stimulus


                            decode();

                            %dp is P(rc(s,t)) = 1 / (num_stimulus * num_transforms)
                            %P(s,s')= P(s'|rc(s,t)) * P(rc(s,t))

                            Psq(s,:) = Psq(s,:)+dp*reshape(pq_r(1:num_stimulus),1,num_stimulus);%probability table
                            Pq(:) = Pq(:) + dp * pq_r(:);%mean assigned probability
                            Qsq(s,best_s) = Qsq(s,best_s)+dp;%frequency table P(s,s')
                            Qq(best_s) = Qq(best_s) + dp;%frequency of each predicted s P(s')
                            %if(s==6)
                            %    Psq
                            %end
                        end
                    end

                    testAvgP = testAvgP + Psq/niter;
                    testAvgQ = testAvgQ + Qsq/niter;

                    %extract info values from frequencies and probabilities

                    for s1 = 1:num_stimulus
                        nb = 0;
                        for s2 = 1:num_stimulus
                            q1 = Qsq(s1,s2);
                            p1 = Psq(s1,s2);

                            %to calculate I with the best matching (not smoothed)
                            if (q1 > eps)
                                Iss2_Q = Iss2_Q + q1 * log2(q1 / (Qq(s2) * Ps));
                                nb = nb + 1;
                            end

                            %to calculate I by probability (smoothed)
                            if (p1 > eps)
                                Iss2_P = Iss2_P + p1 * log2(p1 / (Pq(s2) * Ps));
                            end
                        end
                    end
                    Iss2_Q_avg(nc) =Iss2_Q_avg(nc)+ Iss2_Q / niter;
                    Iss2_P_avg(nc) =Iss2_P_avg(nc)+ Iss2_P / niter;
                end
                testAvgP
                testAvgQ
            end


            MultiCellPlotP(dataIndex,:) = [0 Iss2_P_avg(1:nc_max)];
            MultiCellPlotQ(dataIndex,:) = [0 Iss2_Q_avg(1:nc_max)];
        end


        %uicontrol('Style','text','Position',[100 5 200 20],'String',['num_samples per stim: ' num2str(num_samples*sampleMulti)])

        disp('DONE');

        if (dataIndex==1)
            if(multi_cell_analysis == 1)
                subplot(2,2,1);
                p1= plot([0:1:nc_max],MultiCellPlotP(dataIndex,1:nc_max+1),'k--');
                axis([0 nc_max -0.1 log2(num_stimulus)+0.1]);
                set(p1, 'LineWidth', 1.5);
                hold on

                subplot(2,2,2);
                p2 = plot([0:1:nc_max],MultiCellPlotQ(dataIndex,1:nc_max+1),'k--');
                axis([0 nc_max -0.1 log2(num_stimulus)+0.1]);
                set(p2, 'LineWidth', 1.5);
                hold on
                
                subplot(2,2,3);
            end
            p3 = plot(1:num_cells*num_cells, SingleCellPlot(dataIndex,:),'k--');
            axis([0 num_cells*num_cells -0.1 log2(num_stimulus)+0.1]);
            set(p3, 'LineWidth', 1.5);
            hold on
        elseif (dataIndex==2)
            if(multi_cell_analysis == 1)
                subplot(2,2,1);
                p1= plot([0:1:nc_max],MultiCellPlotP(dataIndex,1:nc_max+1),'k-');
                axis([0 nc_max -0.1 log2(num_stimulus)+0.1]);
                set(p1, 'LineWidth', 1.5);
                hold on

                subplot(2,2,2);
                p2 = plot([0:1:nc_max],MultiCellPlotQ(dataIndex,1:nc_max+1),'k-');
                axis([0 nc_max -0.1 log2(num_stimulus)+0.1]);
                set(p2, 'LineWidth', 1.5);
                hold on

                subplot(2,2,3);
            end
            p3 = plot(1:num_cells*num_cells, SingleCellPlot(dataIndex,:),'k-');
            axis([0 num_cells*num_cells -0.1 log2(num_stimulus)+0.1]);
            set(p3, 'LineWidth', 1.5);
            hold on
        else
            if(multi_cell_analysis == 1)
                subplot(2,2,1);
                p1= plot([0:1:nc_max],MultiCellPlotP(dataIndex,1:nc_max+1),'k:');
                axis([0 nc_max -0.1 log2(num_stimulus)+0.1]);
                set(p1, 'LineWidth', 1.5);
                hold on

                subplot(2,2,2);
                p2 = plot([0:1:nc_max],MultiCellPlotQ(dataIndex,1:nc_max+1),'k:');
                axis([0 nc_max -0.1 log2(num_stimulus)+0.1]);
                set(p2, 'LineWidth', 1.5);
                hold on

                subplot(2,2,3);
            end
            p3 = plot(1:num_cells*num_cells, SingleCellPlot(dataIndex,:),'k:');
            axis([0 num_cells*num_cells -0.1 log2(num_stimulus)+0.1]);
            set(p3, 'LineWidth', 1.5);
            hold on
        end
    end






    function decode()
        eps = 0.0001/num_stimulus;%how to choose the value?
        %eps = 0;
        p_tot =0;
        p_max =0;
        best_s = 0;
        pq_r(1:num_stimulus) = Ps;%start by writing P(s2)
        best_s_Set = [];
        
        
        if (UseMeanForMulti == 0)
            for s2=1:num_stimulus %for each s2, estimate P(r|s2)P(s2)
                
                for c = 1:nc;
                    if(sd(c,s2)<=eps)
                        if(rc(c,s)==ra(c,s2))
                            fact = 1;
                        else
                            fact = eps;
                        end
                    else
                        fact = normpdf(rc(c,s),ra(c,s2),sd(c,s2));
                        if(fact<=eps)
                            fact = eps;
                        end
                    end
                    
                    pq_r(s2) = pq_r(s2) * fact;
                    
                end
                p_tot = p_tot+ pq_r(s2);
                
                %                noise_fact = 1. + eps * (rand() - 0.5);%to randomly choose one from multiple candidates
                if (p_max < pq_r(s2))% * noise_fact)
                    p_max = pq_r(s2);% * noise_fact;
                    %best_s = s2;
                end
            end
        else
            for s2=1:num_stimulus %for each s2, estimate P(r|s2)P(s2)
                for c = 1:nc;
                    
                    
                    
                    
                    
                    if(sd(c,s2)<=eps)
                        if(ra(c,s)==ra(c,s2) && (sd(c,s)<=eps))
                            if(sd(c,s)<=eps)
                                fact = 1;
                            else%katteni tsuketa keshite OK
                                fact = 0;%1/(1+sd(c,s));
                            end%katteni tsuketa keshite OK
                        else
                            fact = eps;
                        end
                    else
                        fact = normpdf(ra(c,s),ra(c,s2),sd(c,s2));
                        if(fact<=eps)%small number threshold
                            fact = eps;
                        end
                    end
                    
                    %test
                    %if(s==6 && s2==6)
                    %    fact
                    %    ra(c,:)
                    %end
                    %test
                    
                    pq_r(s2) = pq_r(s2) * fact;
                    
                end
                p_tot = p_tot+ pq_r(s2);
                
                %                noise_fact = 1. + eps * (rand() - 0.5);%to randomly choose one from multiple candidates
                if (p_max < pq_r(s2))% * noise_fact)
                    p_max = pq_r(s2);% * noise_fact;
                    %best_s = s2;
                end
            end
        end
        
        
        for s2=1:num_stimulus
            if (pq_r(s2) == p_max)
                best_s_Set = [best_s_Set s2];
            end
        end
        
        best_s = best_s_Set(randperm(length(best_s_Set),1));
        
        if p_tot <= eps
            pq_r(1:num_stimulus) = Ps; %if they all died, equal honor to all
            best_s = ceil(num_stimulus * rand());
        else
            pq_r(1:num_stimulus) = pq_r(1:num_stimulus)/ p_tot; % finally normalize to get P(s|r)  */
        end
    end



end

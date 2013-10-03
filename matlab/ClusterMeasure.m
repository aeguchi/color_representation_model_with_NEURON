FR = load('L23_FR_0.txt');
num_stims = length(FR(:,1))
dim = sqrt(length(FR))

localSize = 3;


nTrainings = 21;
c_changes= zeros(nTrainings,1);

listIndex = 1;

for itr = 0:100:2000
    disp(itr)
    FR = load(['L23_FR_', num2str(itr), '.txt']);
    
    
	clusterCoefTot = 0;
	trialCount = 0;
	for stim = 1:num_stims
	    for index = 1:dim*dim
            tmpFR = zeros(localSize*localSize,1);
            for y = 1:localSize-1
                for x = 1:localSize
                    tmpIndex = (y-1)*localSize + x;
                    tmpFR(tmpIndex) = FR(stim,mod(index + (y-1)*dim + (x-1) -1,dim*dim)+1);
                end
            end

            total = 0;
            count = 0;
            for sourceIndex = 1:localSize*localSize
                if sourceIndex==localSize*localSize
                    break
                end
                for targetIndex = sourceIndex+1:localSize*localSize
                    total = total + tmpFR(sourceIndex)*tmpFR(targetIndex);
                    count = count+1;
                end
            end
            clusterCoef = total/count;
            clusterCoefTot = clusterCoefTot + clusterCoef;
            trialCount = trialCount+1;

	    end
	end

	clusterCoefAvg = clusterCoefTot/trialCount
	c_changes(listIndex) = clusterCoefAvg;
    listIndex = listIndex+1;
end
%figure
%hold on
plot([0:100:2000],c_changes, '--ok')
xlabel('training session (300ms/itr)')
ylabel('clustering coefficient')

%FR0 0.3318
%FR100 9.8846

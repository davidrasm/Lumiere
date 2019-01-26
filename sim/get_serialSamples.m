function [sampleTimes, samplesAtTimes, sampleCount align] = get_serialSamples(params, iList, iSeqs)
%   Sample individuals to include in a simulated phylogeny
%   This was formerly get_serialSamples_rhoSampling

birthTimes = iList(:,1);
deathTimes = iList(:,2);
parents = iList(:,3);
states = iList(:,4);
mutations = iList(:,5);
removals = iList(:,6);

align.Names = [];
align.Genotypes = [];

%Get samples
sampleCount = 0;
if (not(params.serialSampling))

    extentLines = find(deathTimes == params.finalTime);
    extentCount = length(extentLines);
    if (params.sampleCount > extentCount)
       display('Sample count to high')
       params.sampleCount = extentCount;
    end

    locs = randsample(extentCount, params.sampleCount);
    linesInSample = num2cell(extentLines(locs));
    linesInSampleArray = cell2mat(linesInSample);

    sampleTimes = params.finalTime;
    samplesAtTimes = cell(0);
    samplesAtTimes{1} = linesInSampleArray;

else
    
        % Serial sampling by node death times
        sampleTimes = [];
        samplesAtTimes = cell(0);
        sortedDeathTimes = sort(unique(deathTimes));
        sampleIndexes = 1:length(iList(:,1));
        for i = 1:length(sortedDeathTimes)
            
            timeSet = sampleIndexes(deathTimes == sortedDeathTimes(i));
            removalSet = sampleIndexes(removals == true);
            timeSet = intersect(timeSet, removalSet);
            totalCount = 0;
            sampleIDs = [];
            
            for m = 1:params.subPops
                
                if (sortedDeathTimes(i) > params.samplingStartTimes(m))
                    
                    popSet = sampleIndexes(states == m);
                    subSet = intersect(popSet,timeSet);
                    
                    if (sortedDeathTimes(i) == params.finalTime)
                        subCount = binornd(length(subSet),params.rho(m));
                    else
                        subCount = binornd(length(subSet),params.samplingFraction(m));
                    end
                    sampleLocs = randsample(length(subSet), subCount);
                    totalCount = totalCount + length(sampleLocs);
                    sampleIDs = [sampleIDs, subSet(sampleLocs)];
                    
                end
                
            end
            
            if (totalCount > 0)
                
                sampleTimes(end+1) = sortedDeathTimes(i);
                samplesAtTimes{end+1} = sampleIDs;
                sampleCount = sampleCount + totalCount;
                
                align.Names = [align.Names; sampleIDs'];
                align.Genotypes = [align.Genotypes; iSeqs(sampleIDs,:)];
                %align.Genotypes = [align.Genotypes; mutations(sampleIDs)];
                
            end
                
        end
    
end

% Somehow lost this sorting step in earlier versions
[sampleTimes,sortedIndexes] = sort(sampleTimes,'descend');
samplesAtTimes = samplesAtTimes(sortedIndexes);


end


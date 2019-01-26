function [times, X, Y, iList, iSeqs, samples, align, simFlag] = run_stochTreeSim_seqEvo(params)

% Simulate birth-death-mutation-sampling process with fitness that depends on seqs/genotypes
% Can also simulate SIR/SIS type dynamics
% Fitness effects can be site-specific as defined by the matrix params.fitMtrx
% For now seqs are represented as a bit vector
% Updated 4.30.18 to allow for sampling to be either coupled or uncoupled from removal

% Params to add
siteList = 1:params.seqLength;
params.initMutClass = 0; % should be zero

simFlag = 0; % goes to one if pop goes extinct

times = 0.0:params.dt:params.finalTime;
totalTimes = length(times);

M = params.subPops;

% Set initial conditions
N = params.popN;
S = params.popN; % should be a column vector
I = zeros(M,1); % should be a column vector
C = zeros(M,M); % cumulative incidence in each pop

if isfield(params, 'initPop')
    initPop = params.initPop;
else
    initPop = randi(M,1); % Randomly seed pop in one subPop
end
I(initPop) = 1;
%S(initPop) = S(initPop) - 1; % hold S constant
        
%Initialize infection tree logger
currInfections{1} = 1;
currInfectionStates{1} = initPop;
iBirthTimes(1) = 0.0;
iDeathTimes(1) = NaN;
iParents(1) = -1;
iStates(1) = initPop;
iMutations(1) = params.initMutClass;
iRemovals(1) = NaN;

% Initialize a new data structure for sequences
% wt genotype is always zero; mut genotype is one
% currMutLoads will just be sum of seq vector
iSeqs = zeros(1,params.seqLength);
%iMutations(1) = sum(iSeqs(1,:)); % in theory

% Initialize currBetas
currMutLoads = sum(iSeqs(1,:)); %params.initMutClass;

%currBetas = params.beta(initPop,:) .* ((1 - params.sd)^currMutLoads(1)); % multiplicative fitness
fitEffects = params.fitMtrx(sub2ind(size(params.fitMtrx), iSeqs+1, siteList));
currBetas = params.beta(initPop,:) * prod(fitEffects); 

% Need currMutRate (one way)
forwardRate = (params.seqLength - currMutLoads(1)) * params.sigma_up;
backwardRate = currMutLoads(1) * params.sigma_down;
currMutRates = forwardRate + backwardRate;

% Added to sample while simulating
sampleTimes = [];
samplesAtTimes = cell(0);
sampleCount = 0;

% Do we need these if we have iSeqs?
align.Names = [];
align.Genotypes = [];

%Run simulation
X = zeros(M,totalTimes,3); % S,I,C
X(:,1,1) = S;
X(:,1,2) = I;

%X(:,1,3) = C;
Y = zeros(M,M,totalTimes);

timeIndex = 2;
nextTime = times(timeIndex);

time = 0.0;

if isfield(params, 'checkPointTime')
    checkPointReached = false;
else
    checkPointReached = true;
end

if isfield(params, 'regimeEndTimes')
    regimeList = 1:1:length(params.regimeEndTimes);
end

while (time < params.finalTime)
    
    % Transmission
    %transRates = params.beta .* (I * (S./N)');
    %transRates(isnan(transRates)) = 0;
    %transRates(isinf(transRates)) = 0;
    %rates(1) = sum(sum(transRates));
    
    totalI = length(currInfections);
    transRates = currBetas .* repmat((S./N)', totalI, 1);
    if (params.logisticGrowth)
        transScaler = 1 - (params.logisticK * params.initSize * exp(params.logisticR * time) / (params.logisticK + params.initSize * (exp(params.logisticR * time)-1)));
        transRates = transScaler * transRates;
    end
    if isfield(params, 'regimeEndTimes')
        future = regimeList(params.regimeEndTimes > time);
        regime = future(1);
        transRates = params.transScaler(regime) * transRates;
    end
    rates(1) = sum(sum(transRates));
    
    % Removal
    recRates = params.nu .* I;
    rates(2) = sum(recRates);
    
    % Host Births
    birthRates = params.mu .* N;
    rates(3) = sum(birthRates); % total M & F birth-rate
    
    % Host Deaths
    deathRates = params.mu .* (I + S);
    rates(4) = sum(deathRates);
    
    % Migrations 
    migRates = params.gamma .* repmat(I,1,M);
    rates(5) = sum(sum(migRates));
    
    % Mutations -- occur independently of population
    %mutRates = (params.sigma_up + params.sigma_down) * sum(I);
    mutRates = sum(currMutRates);
    rates(6) = mutRates;
    
    % Sampling -- can be couple or uncoupled from removal
    if (params.samplingUncoupled)
        samplingRates = params.samplingFraction .* I;
        rates(7) = sum(samplingRates);
    end
    
    [dt, type] = get_nextEventTime(rates);
    time = time + dt;
    
    if (time >= nextTime)
        while (nextTime < time)
            
            %timeIndex
            %I
            
            X(:,timeIndex,1) = S;
            X(:,timeIndex,2) = I;
            
            Y(:,:,timeIndex) = C;
            
            timeIndex = timeIndex + 1;
            if(timeIndex <= totalTimes)
                nextTime = times(timeIndex);
            else
                break
            end
        end
    end
    
    if (time > params.finalTime)
        
        iRemovals(isnan(iDeathTimes)) = true; % so can be sampled
        iDeathTimes(isnan(iDeathTimes)) = params.finalTime;
        
        % Add sampling event
        sampleTimes(end+1) = params.finalTime;
        currI = cell2mat(currInfections);
        currIStates = cell2mat(currInfectionStates);
        sampleIDs = [];
        for m = 1:params.subPops
            popSet = currI(currIStates == m);
            count = binornd(length(popSet),params.rho(m));
            sampleLocs = randsample(length(popSet), count);
            sampleIDs = [sampleIDs, popSet(sampleLocs)];
        end
        samplesAtTimes{end+1} = sampleIDs;
        sampleCount = sampleCount + length(sampleIDs);
        align.Names = [align.Names; sampleIDs'];
        align.Genotypes = [align.Genotypes; iSeqs(sampleIDs,:)];
        
        break;
        
    end
    
    switch(type)
        case 1
            
            % TRANS RATES ARE NOW BY INDIVIDUAL
            parentIndexInCurrInfections = get_nextEvent(sum(transRates,2)); % sums is across columns
            parentIndexInLog = currInfections{parentIndexInCurrInfections};
            parentPop = currInfectionStates{parentIndexInCurrInfections}; % sums is across columns
            childPop = get_nextEvent(transRates(parentIndexInCurrInfections,:));
            
            % Transmission event
            %parentPop = get_nextEvent(sum(transRates,2)); % sums is across columns
            %childPop = get_nextEvent(transRates(parentPop,:));
            %currI = cell2mat(currInfections);
            %currIStates = cell2mat(currInfectionStates);
            %potentialParents = currI(currIStates == parentPop);
            %parentIndexInLog = randsample(potentialParents,1);
            %parentIndexInCurrInfections = find(currI == parentIndexInLog);
            
            %Add infection to infection tree log
            childMutLoad = currMutLoads(parentIndexInCurrInfections);
            childSeq = iSeqs(parentIndexInLog,:);
            childMutRate = currMutRates(parentIndexInCurrInfections);
            add_Infection(time, parentIndexInLog, childPop, childMutLoad, childSeq, childMutRate)
            
            %Update state variables
            %S(childPop) = S(childPop) - 1;
            I(childPop) = I(childPop) + 1;
            %C(childPop) = C(childPop) + 1;
            C(parentPop,childPop) = C(parentPop,childPop) + 1;
            
        case 2
            
            deathPop = get_nextEvent(recRates);
            
            currI = cell2mat(currInfections);
            currIStates = cell2mat(currInfectionStates);
            potentialDeaths = currI(currIStates == deathPop);
            if (length(potentialDeaths) > 1)
                deathIndexInLog = randsample(potentialDeaths,1);
            else
                deathIndexInLog = potentialDeaths(1);
            end
            
            deathIndexInCurrInfections = find(currI == deathIndexInLog); %ceil(rand * I);
            %deathIndexInLog = currInfections{deathIndexInCurrInfections};
            
            % Sample lineage at removal
            if (not(params.samplingUncoupled) && rand < params.samplingFraction(deathPop))
                sampleTimes(end+1) = time;
                samplesAtTimes{end+1} = deathIndexInLog;
                sampleCount = sampleCount + 1;
                align.Names = [align.Names; deathIndexInLog];
                align.Genotypes = [align.Genotypes; iSeqs(deathIndexInLog,:)];
            end
            
            % Remove infection
            removal = true;
            remove_Infection(time, deathIndexInLog, deathIndexInCurrInfections, removal)
            
            %Update state variables
            %if (params.noImmunity(deathPop))
            %    S(deathPop) = S(deathPop) + 1;
            %else
            %    N(deathPop) = N(deathPop) - 1;
            %end
            I(deathPop) = I(deathPop) - 1;
            
        case 3 % Birth (in suceptible pop)
            
            % NOT IMPLEMENTED IN BIRTH-DEATH MODELS
            
        case 4 % Death
            
            % NOT IMPLEMENTED IN BIRTH-DEATH MODELS
            
        case 5 % Migration
            
            parentPop = get_nextEvent(sum(migRates,2)); % sums is over columns
            childPop = get_nextEvent(migRates(parentPop,:));
                
            % Remove random infected from pop
            currI = cell2mat(currInfections);
            currIStates = cell2mat(currInfectionStates);
            potentialTrans = currI(currIStates == parentPop);
            if (length(potentialTrans) > 1)
                transIndexInLog = randsample(potentialTrans,1);
            else
                transIndexInLog = potentialTrans(1);
            end
            transIndexInCurrInfections = find(currI == transIndexInLog); %ceil(rand * I);
            
            childMutLoad = currMutLoads(transIndexInCurrInfections);
            childSeq = iSeqs(transIndexInLog,:);
            childMutRate = currMutRates(transIndexInCurrInfections);
            
            % Remove infection
            removal = false;
            remove_Infection(time, transIndexInLog, transIndexInCurrInfections, removal)
                
            %Add infection to infection tree log
            add_Infection(time, transIndexInLog, childPop, childMutLoad, childSeq, childMutRate)

            %Update state variables
            I(parentPop) = I(parentPop) - 1;
            I(childPop) = I(childPop) + 1;
            
        case 6 % Mutation
            
            % New way with actual seqs
            mutIndexInCurrInfections = get_nextEvent(currMutRates); % had sum(transRates,2) here instead of mutRates
            mutIndexInLog = currInfections{mutIndexInCurrInfections};
            
            % Update seq vec
            currSeq = iSeqs(mutIndexInLog,:);
            siteRates = currSeq * params.sigma_down;
            invSeq = abs(currSeq - 1); % zeros and ones reversed
            siteRates = siteRates + (invSeq * params.sigma_up);
            mutLoc = get_nextEvent(siteRates);
            newSeq = currSeq;
            if (currSeq(mutLoc) == 1)
                newSeq(mutLoc) = 0;
                mutIncrement = -1;
            else
                newSeq(mutLoc) = 1;
                mutIncrement = 1;
            end
            
            % Child here just refers to the new entry in currInfections
            childPop = currInfectionStates{mutIndexInCurrInfections};
            childMutLoad = currMutLoads(mutIndexInCurrInfections) + mutIncrement;
            if (childMutLoad < 0)
                childMutLoad = 0; % make sure mutLoad stays positive
            end
            
            % Compute new mut rates
            forwardRate = (params.seqLength - childMutLoad) * params.sigma_up;
            backwardRate = childMutLoad * params.sigma_down;
            newMutRate = forwardRate + backwardRate;
            
            % Remove infection
            removal = false;
            remove_Infection(time, mutIndexInLog, mutIndexInCurrInfections, removal)
            
            %Add infection to infection tree log
            add_Infection(time, mutIndexInLog, childPop, childMutLoad, newSeq, newMutRate)
            
        case 7 % Sampling
            
            samplePop = get_nextEvent(samplingRates);
            
            currI = cell2mat(currInfections);
            currIStates = cell2mat(currInfectionStates);
            potentialSamples = currI(currIStates == samplePop);
            if (length(potentialSamples) > 1)
                sampleIndexInLog = randsample(potentialSamples,1);
            else
                sampleIndexInLog = potentialSamples(1);
            end
            sampleIndexInCurrInfections = find(currI == sampleIndexInLog); %ceil(rand * I);
            
            % Add sampling event
            sampleTimes(end+1) = time;
            samplesAtTimes{end+1} = sampleIndexInLog;
            sampleCount = sampleCount + 1;
            align.Names = [align.Names; sampleIndexInLog];
            align.Genotypes = [align.Genotypes; iSeqs(sampleIndexInLog,:)];
            
            % Remove infection
            removal = true;
            remove_Infection(time, sampleIndexInLog, sampleIndexInCurrInfections, removal)
            I(samplePop) = I(samplePop) - 1;
            
            
        otherwise
            error('Event type not possible!')    
    end
            
    if (sum(I) <= 0)
        display('Epidemic died off')
        simFlag = 1;
        break;
    end
    
    if (sum(I) > params.maxPopSize)
        display('Pop size over max limit')
        simFlag = 1;
        break;
    end
    
    %if (~checkPointReached) % if not yet reached
    %    if (time >= params.checkPointTime)
    %        display('Passed check point')
    %        if (I < params.checkPointCutoff)
    %            simFlag = 1;
    %            break;
    %        end
    %        checkPointReached = true;
    %    end
    %end
    
end

if (I < 2)
    simFlag = 1;
end

iList = [iBirthTimes', iDeathTimes', iParents', iStates', iMutations', iRemovals'];

[sampleTimes,sortedIndexes] = sort(sampleTimes,'descend');
samplesAtTimes = samplesAtTimes(sortedIndexes);
samples.sampleTimes = sampleTimes;
samples.samplesAtTimes = samplesAtTimes;
samples.sampleCount = sampleCount;


function [dt, type] = get_nextEventTime(rates)

    totalRate = sum(rates);
    dt = log(1/rand) / totalRate;
    eventProbs = cumsum(rates)/totalRate;
    locs = find(eventProbs >= rand);
    type = locs(1);

end

function [index] = get_nextEvent(rates)

    probs = cumsum(rates)/sum(rates);
    locs = find(probs >= rand);
    index = locs(1);

end

function [] = add_Infection(time, parentIndexInLog, childPop, childMutLoad, childSeq, childMutRate)

    %Add infection to tree log
    newIndexInLog = length(iBirthTimes) + 1;
    iBirthTimes(newIndexInLog) = time;
    iDeathTimes(newIndexInLog) = NaN;
    iParents(newIndexInLog) = parentIndexInLog;
    iStates(newIndexInLog) = childPop;
    iMutations(newIndexInLog) = childMutLoad;
    iRemovals(newIndexInLog) = NaN;
    iSeqs(newIndexInLog,:) = childSeq;
    
    % Add infection to currInfections
    newIndex = length(currInfections)+1;
    currInfections{newIndex} = newIndexInLog;
    currInfectionStates{newIndex} = childPop;
    
    % If all mutations carry the same fitness effects
    %currBetas(newIndex,:) = params.beta(childPop,:) .* ((1 - params.sd)^childMutLoad); % multiplicative fitness
    
    % If mutations carry site-specific fitness effects
    fitEffects = params.fitMtrx(sub2ind(size(params.fitMtrx), childSeq+1, siteList));
    currBetas(newIndex,:) = params.beta(childPop,:) * prod(fitEffects);
    
    currMutLoads(newIndex) = childMutLoad;
    currMutRates(newIndex) = childMutRate;
 
end


 
function [] = remove_Infection(time, indexInLog, indexInCurrInfections, removal)

    iDeathTimes(indexInLog) = time;
    iRemovals(indexInLog) = removal;
    
    currInfections(indexInCurrInfections) = [];
    currInfectionStates(indexInCurrInfections) = [];
    currBetas(indexInCurrInfections,:) = [];
    currMutLoads(indexInCurrInfections) = [];
    currMutRates(indexInCurrInfections) = []; % added this for seq evo

end

end





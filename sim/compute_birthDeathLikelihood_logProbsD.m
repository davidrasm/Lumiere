function [L,loggedProbsD,loggedTimes] = compute_birthDeathLikelihood_logProbD(params, phy, genoMap)
% Computes treelikelihood under classic multi-type birth death model
% Tracks each genotype as a different type
% Assumes two-site model for now b/c of mapping between seqs and genotypes

%storeProbs = false;
%trueTipStates = params.trueTipStates;
states = genoMap.genotypes; % all possible genotypes
times = 0:params.dt:params.finalTime;
totalTimes = length(times);
timeIndexes = 1:1:totalTimes;
rescalingCount = 0; % for rescaling probs to avoid numerical underflow
L = -Inf; % set in case function returns at NaN or -Inf prob

% For tracking line probs
loggedProbsD = cell(0);
loggedTimes = cell(0);
    
% Set up other params that should be invariant to genotype state
params.nu = ones(states,1) * params.nu(1);
params.samplingFraction = ones(states,1) * params.samplingFraction(1);
params.rho = ones(states,1) * params.rho(1);

params.beta = diag(genoMap.fitVals * params.beta_0); % only on diagonal
[genoMap,params.gamma] = getTransitionRates(genoMap,params);

% Get envScalers if using logistic growth model
if (params.logisticGrowth)
   envScalers = 1 - (params.logisticK * params.initSize * exp(params.logisticR * times) ./ (params.logisticK + params.initSize * (exp(params.logisticR * times)-1)));
else
   envScalers = ones(1,length(times));
end
    
% First solve probsE backwards in time for all integration times (t=0:params.dt:params.finalTime)
pEInit = ones(states,1) - params.rho;
pEMatrix = solve_probEBack();

%save('mtbd_pEmatrx','pEMatrix');
    
% For each external lineage solve, get tipState and solve lineProbsD
lineProbsD = cell(0);
for lin = 1:phy.tipCount
        
    % Set pD according to sampled state
    currLine = phy.nodes{lin};
    currTime = currLine.lineTimes(1);

    sample_pD = zeros(states,1);
    
    tipSeqVec = currLine.lineSeqs(1,:);
    tipSeq = strcat(num2str(tipSeqVec(1)),num2str(tipSeqVec(2)));
    [genoMap, tipState] = getStateOfSeq(genoMap, tipSeq);
    if (currLine.lineTimes(1) == params.finalTime)
        sample_pD(tipState) = params.rho(tipState);
    else
        sample_pD(tipState) = params.nu(tipState) * params.samplingFraction(tipState); % this should actually be d_i * s_i
    end
    
    if (params.trackLineProbs) % && currLine.save)
       loggedProbsD{lin} = sample_pD; 
       loggedTimes{lin} = currTime;
    end
        
    % Solve for pD along lineage
    lineProbsD{lin} = solve_lineProbsD(currLine, sample_pD);
        
end
    
% For each internal lineage update pD at coalescent event
for lin = phy.tipCount+1:phy.nodeCount
        
    % Get children and their pD's
    currLine = phy.nodes{lin};
    
    if (length(currLine.children) > 1)
        
        % Get init pD for parent
        child1_pD = lineProbsD{currLine.children(1)};
        child2_pD = lineProbsD{currLine.children(2)};

        % This perhaps should be scaledBeta not params.beta
        currTime = currLine.lineTimes(1);
        pastIndexes = timeIndexes(times <= currTime);
        tx = pastIndexes(end); % next time index
        scaledBeta = params.beta * envScalers(tx);

        parent_pD = ((scaledBeta * child1_pD) .* child2_pD) + ((scaledBeta * child2_pD) .* child1_pD);

        if (params.probRescaling)
            while (sum(parent_pD) < params.minThreshProb)
                parent_pD = parent_pD * 10;
                rescalingCount = rescalingCount + 1;
            end
        end

        if (any(sum(isnan(parent_pD))) || any(sum(isinf(parent_pD))))
            %display('Found NaN or infinite updated pD')
            return
        end

        if (params.trackLineProbs) % && currLine.save)
            loggedProbsD{lin} = parent_pD; %[loggedProbsD, parent_pD];
            loggedTimes{lin} = currTime; %[loggedTimes, currTime];
        end
    
    else
        
        parent_pD = lineProbsD{currLine.children(1)};
        
    end
        
    % Solve for pD along lineage
    if (lin == phy.nodeCount)
        lineProbsD{lin} = parent_pD; % this is the root!
    else
        lineProbsD{lin} = solve_lineProbsD(currLine, parent_pD);
    end
        
end
    
% Find root pE
rootTime = phy.nodes{end}.lineTimes(end);
%rootState = phy.nodes{end}.lineStates(end);
pastIndexes = timeIndexes(times <= rootTime);
tx = pastIndexes(end); % next time index

% Compute final likelihood with fixed root state
%root_pE = pEMatrix(rootState,tx);
%root_pD = lineProbsD{end}(rootState);
%L = root_pD / (1 - root_pE);

% Compute final likelihood summing over root state
root_pE = pEMatrix(:,tx);
root_pD = lineProbsD{end};
initFreqs = ones(states,1) / states;
if (params.conditionOnSurvival)
    L = dot(initFreqs, (root_pD ./ (1 - root_pE)));
else
    L = dot(initFreqs, root_pD);
end

if (params.probRescaling)
    logL = log(L) - log(10^rescalingCount);
else
    logL = log(L);
end

L = logL;
    
    function [pEMatrix] = solve_probEBack()
        
        pE = pEInit;
        pEMatrix = zeros(states,totalTimes);
        pEMatrix(:,totalTimes) = pE;
        
        for tx = totalTimes:-1:2
            
            % No time varying parameters
            scaledBeta = params.beta * envScalers(tx);
            
            dtx = times(tx) - times(tx-1);
            
            % death without sampling
            d_pE_1 = (1-params.samplingFraction) .* params.nu;
    
            % no birth or migration
            d_pE_2 = ((scaledBeta + params.gamma) * ones(states,1)) .* pE;
    
            % no death
            d_pE_3 = params.nu .* pE;
    
            % birth of line in j -- both lineages produce no samples
            d_pE_4 = (scaledBeta * pE) .* pE; % simplified
    
            % migration into i from j
            d_pE_5 = params.gamma * pE; % simplified
    
            d_pE = (d_pE_1 - d_pE_2 - d_pE_3 + d_pE_4 + d_pE_5) * dtx;
            pE = pE + d_pE;
            
            pEMatrix(:,tx-1) = pE;
            
        end
        
        if (sum(isnan(pE)) > 0 || sum(isinf(pE)) > 0)
            display('Found negative or infinite pE at time zero')
            keyboard
        end
        
    end

    function [pD] = solve_lineProbsD(currLine, pD)
        
        % Get next line time
        currTime = currLine.lineTimes(1); % sampling time
        endTime = currTime - currLine.parentDistance; % time of birth
        
        % Find next integration time and index tx
        pastIndexes = timeIndexes(times <= currTime);
        tx = pastIndexes(end); % next time index
        nextIntTime = times(tx);
        
        pDMatrix = pD; % store pD values for debugging
            
        while (nextIntTime > endTime) % nextLineTime)

            scaledBeta = params.beta * envScalers(tx);

            pE = pEMatrix(:,tx);

            % no births or migration
            d_pD_1 = ((scaledBeta + params.gamma) * ones(states,1)) .* pD;

            % no deaths
            d_pD_2 = params.nu .* pD;

            % birth of line in j -- line in j produces no samples
            d_pD_3 = (scaledBeta * pE) .* pD;

            % birth of line in j -- line in j produces no samples
            d_pD_4 = (scaledBeta * pD) .* pE;

            % migration into i from j
            d_pD_5 = params.gamma * pD; % simplified 

            d_pD = -(d_pD_1 + d_pD_2) + d_pD_3 + d_pD_4 + d_pD_5;
            dt = times(tx) - times(tx-1);
            pD = pD + (d_pD * dt);
            pDMatrix = [pDMatrix, pD];

            if (params.probRescaling)
                while (sum(pD) < params.minThreshProb)
                    pD = pD * 10;
                    rescalingCount = rescalingCount + 1;
                end
            end

            if (any(sum(isnan(pD))) || any(sum(isinf(pD))))
                %display('Found NaN or infinite updated pD')
                return
            end
            
            if (params.trackLineProbs) % && currLine.save)
                loggedProbsD{lin} = [loggedProbsD{lin}, pD];
                loggedTimes{lin} = [loggedTimes{lin}, times(tx-1)];
            end

            tx = tx - 1;
            nextIntTime = times(tx);

        end
                    
    end

end





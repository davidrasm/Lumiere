function [L,phy] = compute_margFiBirthDeath_trueProbsD_likelihood(params, phy, fitModel,exactLoggedProbsD,loggedTimes)
    % Computes likelihood of phy tree under the marginal fitness
    % birth-death likelihood
    
    % This version uses exact genotype probs stores in exactLoggedProbsD to
    % compute expected fitness of lineages. params.multiplicativeFitness
    % and params.approxEProbs must be set to false for this to work.
    
    % Set up
    states = params.seqStates; % states per site
    sites = params.seqLength; % sites in sequence
    siteList = 1:sites;
    fitClasses = length(fitModel.fitPoints); % fit classes to track in fitSpace 
    if ~isfield(params, 'startTime')
        params.startTime = 0.0;
    end
    times = params.startTime:params.dt:params.finalTime;
    totalTimes = length(times);
    timeIndexes = 1:1:totalTimes;
    L = -Inf; % set return likelihood
    
    % For rescaling probs to avoid numerical underflow
    rescalingCount = 0; 
    fullRescalingCount = 0;
    
    % For tracking line probs
    %loggedProbsD = [];
    %loggedTimes = [];
    
    % Set up fitness effects at each site
    siteFitEffects = params.fitMtrx; % siteFitEffects should come from fitModel
    
    % Set up gamma (mutation rate) matrix
    upVec = (ones(states-1, 1) * params.sigma_up);
    downVec = (ones(states-1, 1) * params.sigma_down);
    params.gamma = diag(upVec, 1) + diag(downVec, -1);
    
    % Set up other params that should be invariant to seq changes
    params.nu = ones(states,1) * params.nu(1);
    params.samplingFraction = ones(states,1) * params.samplingFraction(1);
    params.rho = ones(states,1) * params.rho(1);
    
    % Get envScalers if using logistic growth model
    if (params.logisticGrowth)
       envScalers = 1 - (params.logisticK * params.initSize * exp(params.logisticR * times) ./ (params.logisticK + params.initSize * (exp(params.logisticR * times)-1)));
    else
       envScalers = ones(1,length(times));
    end
    
    % First solve probsE backwards in time for all integration times (t=0:params.dt:params.finalTime)
    pEMatrix = solve_probE_fitSpace();
    
    %mfbd_pEMatrix = pEMatrix;
    %save('mfbd_pEmatrx','mfbd_pEMatrix');
    
    % For each external lineage, get tipState and solve siteLineProbsD
    siteLineProbsD = cell(0); % holds marginal site prob densities - or D_n,k,i
    fullLineProbsD = cell(0); % holds whole lineage prob densities - or D_n
    for lin = 1:phy.tipCount
        
        % Set site pDs according to sampled state
        currLine = phy.nodes{lin};
        sample_pD = zeros(states,sites);
        tipSeq = currLine.lineSeqs(1,:) + 1; % + 1 so indexed from 1
        if (isempty(currLine.lineTimes))
            lineStartTime = currLine.date; % if loading an external phylogeny
        else
            lineStartTime = currLine.lineTimes(1); % if lineTimes are stored from simulation
        end
        if (lineStartTime == params.finalTime)
            for s = 1:sites
                sample_pD(tipSeq(s),s) = params.rho(tipSeq(s)); % sampling probability is rho at final time
            end
            full_sample_pD = params.rho(1);
        else
            for s = 1:sites
                sample_pD(tipSeq(s),s) = params.nu(tipSeq(s)) * params.samplingFraction(tipSeq(s)); % this should actually be d_i * s_i
            end
            full_sample_pD = params.nu(1) * params.samplingFraction(1);
        end
        
        if (params.trackLineProbs)
            [genotypeProbs] = get_genotypeProbs(sample_pD); % or norm_pD
            fitValue = dot(genotypeProbs,fitModel.genotypeFitness); % genotypeFitness
            phy.nodes{lin}=pushLineEvent(phy.nodes{lin},fitValue,zeros(1,sites),lineStartTime);
        end
        
        % Solve for pD along lineage
        [siteLineProbsD{lin}, fullLineProbsD{lin}] = solve_siteLineProbsD(currLine, sample_pD, full_sample_pD);
        
    end
    
    % For each internal lineage update pD at coalescent event
    for lin = phy.tipCount+1:phy.nodeCount
        
        currLine = phy.nodes{lin};
        
        if (length(currLine.children) > 1)
        
            if (isempty(currLine.lineTimes))
                currTime = currLine.date;
            else
                currTime = currLine.lineTimes(1);
            end
            pastIndexes = timeIndexes(times <= currTime);
            tx = pastIndexes(end); % next time index
            scaledBeta = params.beta_0 * envScalers(tx);

            % Get init pD for parent
            child1_pD = siteLineProbsD{currLine.children(1)};
            child2_pD = siteLineProbsD{currLine.children(2)};

            % Get normalized site pD for both children
            norm_child1_pD = child1_pD ./ repmat(sum(child1_pD),states,1); % normalize by dividing each column by its sum
            norm_child2_pD = child2_pD ./ repmat(sum(child2_pD),states,1); % normalize by dividing each column by its sum
            
            if (params.multiplicativeFitness)
            
                weightedEffects_child1 = sum(norm_child1_pD .* siteFitEffects); % fitness effects weighted by their probs
                expFitness = prod(weightedEffects_child1) ./ weightedEffects_child1; % take product over all other sites
                marginalFitEffects_child1 = siteFitEffects .* repmat(expFitness,states,1);

                weightedEffects_child2 = sum(norm_child2_pD .* siteFitEffects); % fitness effects weighted by their probs
                expFitness = prod(weightedEffects_child2) ./ weightedEffects_child2; % take product over all other sites
                marginalFitEffects_child2 = siteFitEffects .* repmat(expFitness,states,1);
                
                child1_full_Beta = scaledBeta * prod(weightedEffects_child1);
                child2_full_Beta = scaledBeta * prod(weightedEffects_child2);
            
            else
                
                loggedGenotypeProbs_child1 = exactLoggedProbsD{currLine.children(1)}(:,end);
                loggedGenotypeProbs_child2 = exactLoggedProbsD{currLine.children(2)}(:,end);

                marginalFitEffects_child1 = zeros(states,sites);
                marginalFitEffects_child2 = zeros(states,sites);
                for site = 1:sites
                    for state = 1:states
                       %[genotypeProbs] = get_condGenotypeProbs(norm_child1_pD,site,state);
                       [genotypeProbs] = get_exactCondGenotypeProbs(loggedGenotypeProbs_child1,site,state);
                       marginalFitEffects_child1(state,site) = dot(genotypeProbs,fitModel.genotypeFitness);
                       %[genotypeProbs] = get_condGenotypeProbs(norm_child2_pD,site,state);
                       [genotypeProbs] = get_exactCondGenotypeProbs(loggedGenotypeProbs_child2,site,state);
                       marginalFitEffects_child2(state,site) = dot(genotypeProbs,fitModel.genotypeFitness);
                    end
                end
                
                %[genotypeProbs] = get_genotypeProbs(norm_child1_pD);
                loggedGenotypeProbs_child1 = loggedGenotypeProbs_child1 ./ sum(loggedGenotypeProbs_child1);
                child1_full_Beta = scaledBeta * dot(loggedGenotypeProbs_child1,fitModel.genotypeFitness);
                %[genotypeProbs] = get_genotypeProbs(norm_child2_pD);
                loggedGenotypeProbs_child2 = loggedGenotypeProbs_child2 ./ sum(loggedGenotypeProbs_child2);
                child2_full_Beta = scaledBeta * dot(loggedGenotypeProbs_child2,fitModel.genotypeFitness);
            
            end
            
            parent_pD = zeros(states,1);        
            for s = 1:sites
               child1_Beta = scaledBeta .* diag(marginalFitEffects_child1(:,s)); % accounts for fitness effects at this site and the expected effect at all other sites
               child2_Beta = scaledBeta .* diag(marginalFitEffects_child2(:,s)); % accounts for fitness effects at this site and the expected effect at all other sites
               parent_pD(:,s) = ((child2_Beta * child1_pD(:,s)) .* child2_pD(:,s)) + ((child1_Beta * child2_pD(:,s)) .* child1_pD(:,s));
            end

            % Update full_pD
            child1_full_pD = fullLineProbsD{currLine.children(1)};
            child2_full_pD = fullLineProbsD{currLine.children(2)};
            parent_full_pD = (child2_full_Beta * child1_full_pD * child2_full_pD) + (child1_full_Beta * child2_full_pD * child1_full_pD);
            
            % Rescale to prevent under/overflow
            if (params.probRescaling)
                [parent_pD, parent_full_pD] = rescale_densities(parent_pD, parent_full_pD);
            end

            if (any(sum(isnan(parent_pD))) || any(sum(isinf(parent_pD))))
                display('Found NaN or infinite updated pD')
                return
            end

            if (params.trackLineProbs)                
                [genotypeProbs] = get_genotypeProbs(parent_pD); % or norm_pD
                fitValue = dot(genotypeProbs,fitModel.genotypeFitness);
                phy.nodes{lin}=pushLineEvent(phy.nodes{lin},fitValue,zeros(1,sites),currTime);
            end
        
        else
            
            parent_pD = siteLineProbsD{currLine.children(1)};
            parent_full_pD = fullLineProbsD{currLine.children(1)};
            
        end
        
        % Solve for pD along lineage
        if (lin == phy.nodeCount)
            siteLineProbsD{lin} = parent_pD; % this is the root!
            fullLineProbsD{lin} = parent_full_pD;
        else
            [siteLineProbsD{lin}, fullLineProbsD{lin}] = solve_siteLineProbsD(currLine, parent_pD, parent_full_pD);
        end
        
    end
    
    % Find root pE
    if isempty(phy.nodes{end}.lineTimes)
        rootTime = phy.nodes{end}.date; % not sure what we should put here
    else
        rootTime = phy.nodes{end}.lineTimes(end);
    end
    pastIndexes = timeIndexes(times <= rootTime);
    tx = pastIndexes(end); % next time index
    
    root_pD = siteLineProbsD{end};
    norm_pD = root_pD ./ repmat(sum(root_pD),states,1); % normalize by dividing each column by its sum
    
    if (params.approxEProbs)
                
        % What if fitness effects are not multiplicative?
        weightedEffects = sum(norm_pD .* siteFitEffects); % fitness effects weighted by their probs
        expFitness = prod(weightedEffects) ./ weightedEffects; % take product over all other sites
        marginalFitEffects = siteFitEffects .* repmat(expFitness,states,1); % site fitness effects marginalized over all other sites
    
        % Compute site-specific pE values from fitness landscape
        root_pE = zeros(states,sites);
        for site = 1:sites
            for state = 1:states
                fitValue = marginalFitEffects(state,site);
                [fitModel,eFitClass] = getDiscreteFitPoint(fitModel,fitValue);
                root_pE(state,site) = pEMatrix(eFitClass,tx);
            end
        end
    
        fitValue = prod(weightedEffects);
        [fitModel,eFitClass] = getDiscreteFitPoint(fitModel,fitValue);
        root_line_pE = pEMatrix(eFitClass,tx);
                
    else
        
        root_pE = zeros(states,sites);
        for site = 1:sites
            for state = 1:states 
                [genotypeProbs] = get_condGenotypeProbs(norm_pD,site,state);
                root_pE(state,site) = dot(genotypeProbs,pEMatrix(:,tx));
            end
        end
        
        [genotypeProbs] = get_genotypeProbs(norm_pD);
        root_line_pE = dot(genotypeProbs,pEMatrix(:,tx));       
                
    end
    
    % Need to get root_full_pD and root_full_pE;    
    if (params.conditionOnSurvival)
        full_tree_pD = fullLineProbsD{end} / (1 - root_line_pE);
    else
        full_tree_pD = fullLineProbsD{end};
    end
    
    % Compute final likelihood
    initFreqs = ones(states,1) / states; % assuming a uniform prior
    logL = 0;
    for s = 1:sites
        if (params.conditionOnSurvival)
            logL = logL + log(dot(initFreqs, (root_pD(:,s) ./ (1 - root_pE(:,s)))));
        else
            logL = logL + log(dot(initFreqs, root_pD(:,s)));
        end
    end
    
    if (params.probRescaling)
        scaled_logL = logL - log(10^rescalingCount);
        scaled_full_tree_pD = log(full_tree_pD) - log(10^fullRescalingCount);
        logL = scaled_logL - ((sites-1) * scaled_full_tree_pD);
    else
        logL = logL - ((sites-1) * log(full_tree_pD)); % 'divide' by full prob of tree
    end
    
    L = logL;
    
    function [pEMatrix] = solve_probE_fitSpace()
    
        % Set up params in fitness space
        params.eNu = ones(fitClasses,1) * params.nu(1);
        params.eSamplingFraction = ones(fitClasses,1) * params.samplingFraction(1);
        params.eRho = ones(fitClasses,1) * params.rho(1);
        params.eBeta = params.beta_0 * diag(fitModel.fitPoints);
        [fitModel, params.eGamma] = getTransitionRates(fitModel, params);
        
        % Set up pEMatrix to store E values
        pEInit = ones(fitClasses,1) - params.eRho;
        pE = pEInit;
        pEMatrix = zeros(fitClasses,totalTimes);
        pEMatrix(:,totalTimes) = pE;
        
        for tx = totalTimes:-1:2
            
            scaledBeta = params.eBeta * envScalers(tx);
            
            dtx = times(tx) - times(tx-1);
            
            % death without sampling
            d_pE_1 = (1-params.eSamplingFraction) .* params.eNu;
    
            % no birth or migration
            d_pE_2 = ((scaledBeta + params.eGamma) * ones(fitClasses,1)) .* pE;
    
            % no death
            d_pE_3 = params.eNu .* pE;
    
            % birth of line in j -- both lineages produce no samples
            d_pE_4 = (scaledBeta * pE) .* pE; % simplified
    
            % migration into i from j
            d_pE_5 = params.eGamma * pE; % simplified
    
            d_pE = (d_pE_1 - d_pE_2 - d_pE_3 + d_pE_4 + d_pE_5) * dtx;
            pE = pE + d_pE;
            
            pEMatrix(:,tx-1) = pE;
            
        end
        
        if (sum(isnan(pE)) > 0 || sum(isinf(pE)) > 0)
            display('Found negative or infinite pE at time zero')
            return
        end
        
    end

    function [pD, full_pD] = solve_siteLineProbsD(currLine, pD, full_pD)
        
        % Get start and end time for lineage
        if (isempty(currLine.lineTimes))
            currTime = currLine.date;
        else
            currTime = currLine.lineTimes(1);
        end
        endTime = currTime - currLine.parentDistance; % time of birth
        
        % Find next integration time and index tx
        pastIndexes = timeIndexes(times <= currTime);
        tx = pastIndexes(end); % next time index
        nextIntTime = times(tx); % next integration time
        
        % Get exact genotype probs
        %loggedGenotypeProbs = exactLoggedProbsD{lin}(:,1); % just for
        %debugging
        %loggedTime = loggedTimes{lin}(1); % just for
        %debugging
        loggedIndex = 1;
            
        while (nextIntTime > endTime)
                
            scaledBeta = params.beta_0 * envScalers(tx);
            
            % Compute normalized marginal site pDs (omega_n,k,i)
            norm_pD = pD ./ repmat(sum(pD),states,1); % normalize by dividing each column by its sum
            
            % Get true (logged) genotype probs
            loggedGenotypeProbs = exactLoggedProbsD{lin}(:,loggedIndex);
            normalizedGenotypeProbs = loggedGenotypeProbs ./ sum(loggedGenotypeProbs);
            loggedIndex = loggedIndex + 1;
            
            if (params.multiplicativeFitness)
               
                % Can compute marginalFitEffects for all sites at once (if no epistasis)
                weightedEffects = sum(norm_pD .* siteFitEffects); % fitness effects weighted by their probs
                expFitness = prod(weightedEffects) ./ weightedEffects; % take product over all other sites
                %if (params.coupleFitnessEffects) % can uncouple fitness effects at other sites
                    marginalFitEffects = siteFitEffects .* repmat(expFitness,states,1); % site fitness effects marginalized over all other sites
                %else
                    %marginalFitEffects(state,site) = siteFitEffects(state,site);
                %end
                expBeta = scaledBeta * prod(weightedEffects); % expected beta for entire lineage
                
            else
                
                % Compute marginalFitEffects based on conditional geno probs
                marginalFitEffects = zeros(states,sites);
                for site = 1:sites
                    for state = 1:states 
                        %[genotypeProbs] = get_condGenotypeProbs(norm_pD,site,state);
                        [genotypeProbs] = get_exactCondGenotypeProbs(loggedGenotypeProbs,site,state);
                        %if (params.coupleFitnessEffects)
                            marginalFitEffects(state,site) = dot(genotypeProbs,fitModel.genotypeFitness);
                        %else
                            %marginalFitEffects(state,site) = siteFitEffects(state,site);
                        %end
                    end
                end
                %[genotypeProbs] = get_genotypeProbs(norm_pD);
                expBeta = scaledBeta * dot(normalizedGenotypeProbs,fitModel.genotypeFitness);
                
            end
            
            % Get pE
            pE = zeros(states,sites);
            if (params.approxEProbs)
                
                % Compute site-specific pE values based on current fitness
                for site = 1:sites
                    for state = 1:states
                        fitValue = marginalFitEffects(state,site);
                        [fitModel,eFitClass] = getDiscreteFitPoint(fitModel,fitValue);
                        pE(state,site) = pEMatrix(eFitClass,tx);
                    end
                end
                
                % Find pE value for entire lineage
                fitValue = prod(weightedEffects);
                [fitModel,eFitClass] = getDiscreteFitPoint(fitModel,fitValue);
                line_pE = pEMatrix(eFitClass,tx);
                
            else
                
                % Compute pE values based on conditional genotype probs
                for site = 1:sites
                    for state = 1:states
                        %[genotypeProbs] = get_condGenotypeProbs(norm_pD,site,state);
                        [genotypeProbs] = get_exactCondGenotypeProbs(loggedGenotypeProbs,site,state);
                        pE(state,site) = dot(genotypeProbs,pEMatrix(:,tx));
                    end
                end
                
                % For entire lineage
                %[genotypeProbs] = get_genotypeProbs(norm_pD);
                line_pE = dot(normalizedGenotypeProbs,pEMatrix(:,tx));
                
            end
            
            if (params.trackLineProbs)
                [genotypeProbs] = get_genotypeProbs(norm_pD);
                fitValue = dot(genotypeProbs,fitModel.genotypeFitness);
                phy.nodes{lin}=pushLineEvent(phy.nodes{lin},fitValue,zeros(1,sites),times(tx));
            end
            
            % Update pD
                
            % no births or migration
            d_pD_1 = ((scaledBeta * marginalFitEffects) + (params.gamma * ones(states,sites))) .* pD;
                
            % no deaths
            d_pD_2 = repmat(params.nu,1,sites) .* pD;
                
            % birth of line in j -- line in j produces no samples
            d_pD_3 = (scaledBeta * marginalFitEffects) .* pE .* pD;

            % birth of line in j -- line in j produces no samples
            d_pD_4 = (scaledBeta * marginalFitEffects) .* pD .* pE; % same as above?

            % migration into i from j
            d_pD_5 = params.gamma * pD;% without loops
                
            d_pD = -(d_pD_1 + d_pD_2) + d_pD_3 + d_pD_4 + d_pD_5;
            dt = times(tx) - times(tx-1);
            pD = pD + (d_pD * dt);
            
            % Update full_pD for entire lineage
            d_full_pD_1 = (expBeta + params.nu(1)) * full_pD;
            d_full_pD_2 = 2 * expBeta * line_pE * full_pD; % assumes pE's are all equal
            d_full_pD = -d_full_pD_1 + d_full_pD_2;
            full_pD = full_pD + (d_full_pD * dt);
            
            % Rescale to prevent under/overflow
            if (params.probRescaling)
                [pD, full_pD] = rescale_densities(pD,full_pD);
            end

            if (any(sum(isnan(pD))) || any(sum(isinf(pD))))
                display('Found NaN or infinite updated pD')
                return
            end
            
            tx = tx - 1;
            nextIntTime = times(tx);
            
        end
                    
    end

    function [pD, full_pD] = rescale_densities(pD,full_pD)
        while (any(sum(pD) < params.minThreshProb))
            reScaleSites = siteList(sum(pD) < params.minThreshProb);
            pD(:,reScaleSites) = pD(:,reScaleSites) * 10;
            rescalingCount = rescalingCount + length(reScaleSites);
        end
        while (any(sum(pD) > params.maxThreshProb))
            reScaleSites = siteList(sum(pD) > params.maxThreshProb);
            pD(:,reScaleSites) = pD(:,reScaleSites) / 10;
            rescalingCount = rescalingCount - length(reScaleSites);
        end
        while (full_pD < params.minThreshProb)
            full_pD = full_pD * 10;
            fullRescalingCount = fullRescalingCount + 1;
        end
        while (any(sum(full_pD) > params.maxThreshProb))
            full_pD = full_pD / 10;
            fullRescalingCount = fullRescalingCount - 1;
        end 
    end

    function [genotypeProbs] = get_genotypeProbs(norm_pD)
       
        genotypeProbs = zeros(fitModel.genotypes,1);
        for g = 1:fitModel.genotypes
            margSiteProbs = zeros(1,sites);
            for s = 1:sites
                %type = str2num(genoMap.seqArray{g}(s)) + 1;
                type = fitModel.genoSeqArray(g,s) + 1;
                margSiteProbs(s) = norm_pD(type,s);
            end
            genotypeProbs(g) = prod(margSiteProbs); % will these be normalized? no   
        end
        if (sum(genotypeProbs) <= 0)
            genotypeProbs = ones(fitModel.genotypes,1) / fitModel.genotypes;
        else
            genotypeProbs = genotypeProbs / sum(genotypeProbs); % Renormalize
        end
        
    end

    function [genotypeProbs] = get_condGenotypeProbs(norm_pD,site,state)
       
        % Compute genotype probs conditional on knowing state at one site
        
        genotypeProbs = zeros(fitModel.genotypes,1);
        for g = 1:fitModel.genotypes
            %if (str2num(genoMap.seqArray{g}(site)) + 1 ~= state)
            if (fitModel.genoSeqArray(g,site) + 1 ~= state)
                genotypeProbs(g) = 0; 
            else
                margSiteProbs = zeros(1,sites);
                for s = 1:sites
                    if s ~= site
                        %type = str2num(genoMap.seqArray{g}(s)) + 1;
                        type = fitModel.genoSeqArray(g,s) + 1;
                        margSiteProbs(s) = norm_pD(type,s);
                    else
                        margSiteProbs(s) = 1;
                    end
                end
                genotypeProbs(g) = prod(margSiteProbs); % will these be normalized? no
            end
        end
        if (sum(genotypeProbs) <= 0)
            genotypeProbs = ones(fitModel.genotypes,1) / fitModel.genotypes;
        else
            genotypeProbs = genotypeProbs / sum(genotypeProbs); % Renormalize
        end
        
    end

    function [genotypeProbs] = get_exactCondGenotypeProbs(genotypeProbs,site,state)
       
        % Compute genotype probs conditional on knowing exact genotype
        % probs
        
        %genotypeProbs = zeros(fitModel.genotypes,1);
        %genotypeProbs = genotypeProbs ./ sum(genotypeProbs);
        for g = 1:fitModel.genotypes
            if (fitModel.genoSeqArray(g,site) + 1 ~= state)
                genotypeProbs(g) = 0; 
            end
        end
        if (sum(genotypeProbs) <= 0)
            genotypeProbs = ones(fitModel.genotypes,1) / fitModel.genotypes;
        else
            genotypeProbs = genotypeProbs / sum(genotypeProbs); % Renormalize
        end
        
    end

    function [siteProb] = get_margSiteProb(genotypeProbs,site,state)
       
        % Compute marginal site probs given joint genotype probs
        siteProb = 0;
        for g = 1:fitModel.genotypes
            %if (str2num(genoMap.seqArray{g}(site)) + 1 == state)
            if (fitModel.genoSeqArray(g,site) + 1 == state)
                siteProb = siteProb + genotypeProbs(g); 
            end
        end
        
    end


end





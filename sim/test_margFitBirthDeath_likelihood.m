function [ output_args ] = test_margFitBirthDeath_likelihood( input_args )
% Test implemenation of marginal fitness birth death model by computing
% likelihood profiles for simulated phylogenies

% Set simulations  control params
params.dt = 0.1; % sim/integration time step
params.finalTime = 100.0; %stop time
params.trueTipStates = true; % assign tips true states at sampling times?
params.displayTrees = true; % display simulated tree?
params.trackLineProbs = true; % track probabilites along lineages
params.minTreeSize = 20; % minimum threshold size for simulated tree
params.maxTreeSize = 60; % maximum threshold size for simulated tree

% Set options for computing BD likelihood
params.approxEProbs = false; % approximate E probs?
params.probRescaling = true; % rescale prob densities to avoid numerical under/overflow? should be true
params.minThreshProb = 1e-200; % for rescaling prob densities
params.maxThreshProb = 1e200; % for rescaling prob densities
params.multiplicativeFitness = false; % if fitness effects are multiplicative, can make likelihood calcs faster
params.conditionOnSurvival = true; % condition likelihood on survival? should be true
params.coupleFitnessEffects = true; % now assumed to always be true

% Mutation/fitness params
params.initMutClass = 0;
params.seqStates = 2; % number of states at each site
params.seqLength = 2; % number of sites
params.sd = 0.5; % selective cost of mutations: fitness = (1-sd)^muts
params.epiFitness = (1-params.sd)^2; % fitness of double mutant with possible epistasis
params.sigma_up = 0.05; % forward mut rate
params.sigma_down = 0.05; % back mut rate

% Set up fitness model
params.fitMtrx = ones(2,params.seqLength);
params.fitMtrx(2,1:params.seqLength) = 1 - params.sd;
fitModel = fitnessModel(params.fitMtrx,params);

% Demographic rates
params.beta_0 = 0.25; % base birth rate (lambda_0 in MFBD paper) 
params.nu = 0.05; % Death/Removal rates (d in MFBD paper)

% Environmental params
params.logisticGrowth = true; % rescales birth rates to get logistic growth
params.logisticK = 0.5;
params.logisticR = 0.5;
params.initSize = 0.005;

% Parameters only used for tree simulations:
params.subPops = 1; % should be states
params.noImmunity = ones(params.subPops,1); % if true, then SIS; if false, then SIR
params.beta = diag(ones(params.subPops,1) * params.beta_0); % birth rates
params.gamma = ones(params.subPops,params.subPops) * 0.0; % migration rates in sim
params.mu = zeros(params.subPops,1); % Host birth/death rates
genoMap = epistaticGenotypeMap(params.fitMtrx,params); % Still need this to simulate phylogenies
params.initPop = 1; 
params.popN = ones(params.subPops,1);

% Sampling params
params.serialSampling = true;
params.samplingUncoupled = false; % has to be false, never implmented uncoupled sampling in run_stochTreeSim_genotypeEvo 
params.samplingStartTimes = ones(params.subPops,1) * 0.01;
params.samplingFraction = ones(params.subPops,1) * 0.05; % sampling in past
params.rho = ones(params.subPops,1) * 0.05; % sampling at present

% Get a tree by forward Gillespie simulation
simFlag = 1;
while(simFlag)
    
    [times, X, Y, iList, iSeqs, simFlag] = run_stochTreeSim_genotypeEvo(params, genoMap); %Run stochastic sim
    [sampleTimes, samplesAtTimes, sampleCount, align] = get_serialSamples(params, iList, iSeqs);
    
    if (sampleCount < params.minTreeSize || sampleCount > params.maxTreeSize)
        strout = strcat('Tree contained', {' '}, num2str(sampleCount), {' '}, 'samples');
        display(strout)
        simFlag = 1;
    end
    
end

% Plot sim pop dynamics
figure();
fTime = params.finalTime / params.dt;
hmap = colormap(parula(4));
for m = 1:params.subPops
    plot(times(1:fTime), X(m,1:fTime,2), '-', 'LineWidth', 2.0, 'Color', hmap(m,:)); hold on;
end
box off; xlabel('Time'); ylabel('Prevalence');

% Construct a new phylo object from iList and the samples
iList(:,5) = iList(:,5) + 1; % shift mutLoad count
phy = phylo(iList,iSeqs,sampleTimes,samplesAtTimes);
phy = convertSeqsToStates(phy,genoMap); % Converts seqs to states in genotypeMap
if (params.displayTrees)
    fileName = 'simPhyTree.tre';
    [phy, treeString] = writeSimMapTree(phy, fileName);
    mut_hmap = colormap(parula(4));
    mutTreeParams.colorMap = mut_hmap;
    mutTreeParams.colorBar = false;
    simmaptree_plot(treeString, mutTreeParams);
end

% Check likelihood under exact MTBD model tracking genotypes
%[trueL,exactLoggedProbsD,loggedTimes] = compute_birthDeathLikelihood_logProbsD(params, phy, genoMap);
%[L,phy] = compute_margFitBirthDeath_trueProbsD_likelihood(params, phy, fitModel, exactLoggedProbsD, loggedTimes);

% Compute entire likelihood surface
trueValY = params.beta_0;
trueValX =  params.sd; %params.epiFitness;
paramVecY = trueValY;
paramVecX = 0.0:0.05:0.95; % sd/epiFitness
paramVecX = [paramVecX, 0.99];

xVals = length(paramVecX);
yVals = length(paramVecY);
likeMatrix = zeros(yVals,xVals);
trueLikeMatrix = zeros(yVals,xVals);
condLikeMatrix = zeros(yVals,xVals);

totalReps = xVals * yVals
for x = 1:xVals
    for y = 1:yVals
        
        rep = (yVals*(x-1)) + y
        
        % Copy new params to paramsCopy
        paramsCopy = params;
        paramsCopy.beta_0 = paramVecY(y);
        
        % If estimating beta
        %paramsCopy.beta_0 = paramVecX(x);
        %paramsCopy.beta = diag(ones(params.subPops,1) * paramsCopy.beta_0); % birth rates
        
        % If estimating sd
        paramsCopy.sd = paramVecX(x);
        paramsCopy.fitMtrx(2,1:params.seqLength) = 1 - paramsCopy.sd;
        paramsCopy.epiFitness = (1-paramsCopy.sd)^2; %1 - paramVecX(x);  
        
        % Update genotypeMap and fitnessModel
        [genoMap] = updateGenotypeMap(genoMap,paramsCopy.fitMtrx, paramsCopy);
        [fitModel] = updateFitLandscape(fitModel, paramsCopy.fitMtrx, paramsCopy);
        
        % Compute exact MTBD likelihood tracking each genotype
        %trueL = compute_birthDeathLikelihood(paramsCopy, phy, genoMap);
        [trueL,exactLoggedProbsD,loggedTimes] = compute_birthDeathLikelihood_logProbsD(paramsCopy, phy, genoMap);
        
        tic
        
        % Compute approx MFBD likelihood
        [L,phy] = compute_margFitBirthDeath_likelihood(paramsCopy, phy, fitModel);
        
        % Conditioned on exact genotype probs
        [condL,phy] = compute_margFitBirthDeath_trueProbsD_likelihood(paramsCopy, phy, fitModel, exactLoggedProbsD, loggedTimes); 
        
        toc
        
        trueLikeMatrix(y,x) = real(trueL);
        likeMatrix(y,x) = real(L);
        condLikeMatrix(y,x) = real(condL);
        

    end
end

save('fitBD_simpleTest_multiFitness', 'params', 'phy', 'paramVecX', 'paramVecY', 'likeMatrix', 'trueLikeMatrix');

figure()
minValue = max(max(likeMatrix)) - 50;
likeMatrix(likeMatrix < minValue) = minValue;
trueLikeMatrix(trueLikeMatrix < minValue) = minValue;
condLikeMatrix(condLikeMatrix < minValue) = minValue;
plot(paramVecX, likeMatrix', 'Color', 'b', 'LineWidth', 2.0); hold on;
plot(paramVecX, trueLikeMatrix', 'Color', 'k', 'LineWidth', 2.0); hold on;
plot(paramVecX, condLikeMatrix', 'Color', 'r', 'LineWidth', 2.0); hold on;
line([trueValX, trueValX], ylim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'b')
xlabel('Epistatic fitness cost \epsilon','FontSize',14); 
ylabel('Log likelihood','FontSize',14);

% If 2D
% figure()
% minValue = max(max(likeMatrix)) - 100;
% likeMatrix(likeMatrix < minValue) = minValue; 
% surf(paramVecX,paramVecY,likeMatrix);
% line([trueValX, trueValX], [trueValY, trueValY], zlim, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'k')

end


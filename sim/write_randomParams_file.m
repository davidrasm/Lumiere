function [ output_args ] = write_randomParams_file(file_out)
% Write params to file for tree simulations
% Random values are generataed for key parameters

if ~exist('file_out','var')
  file_out = 'fitBD_neutralFitModel_5sites';
end

%Set epi params
params.subPops = 1;
params.dt = 0.1; %sim time step
params.finalTime = 100.0; %stop time
params.noImmunity = ones(params.subPops,1); % if true, then SIS; if false, then SIR

% Set control params
params.trueTipStates = true; % assign tips true states at sampling times?
params.trackLineProbs = false;
params.displayTrees = false;
params.minTreeSize = 100;
params.maxTreeSize = 500;
params.maxPopSize = 10000; % used to kill sims that have grown to large

% Sampling params
params.serialSampling = true; % still used anywhere?
params.samplingUncoupled = false;
params.samplingStartTimes = ones(params.subPops,1) * 0.01;
randSamplingRate = 0.1 + (1-0.1)*rand(1,1); % Uniform[0.1 - 1]
params.samplingFraction = ones(params.subPops,1) * randSamplingRate; % also the samplingProportion if sampling is uncoupled
params.rho = ones(params.subPops,1) * randSamplingRate;

% Pop order:
params.initPop = 1;
params.popN = ones(params.subPops,1); % just a place holder

% Mutation params
params.initMutClass = 0;
params.seqLength = 2; % was previously params.maxMutClass
randMutRate = exprnd(0.02); % exponential to get mostly lower values
params.sigma_up = randMutRate; % forward mut rate
params.sigma_down = randMutRate; % back mut rate

params.fitMtrx = ones(2,params.seqLength);

% With constant fitness costs 
%params.sd = 0.25; % selective cost of mutations: fitness = (1-sd)^muts
%params.fitMtrx(2,1:params.seqLength) = 1 - params.sd;

% With random fitness costs (lognormal)
%trueSiteEffects = lognrnd(-0.2231, 0.2, 1, params.seqLength);
%trueSiteEffects(trueSiteEffects > 1.0) = 1.0;

% With random fitness effects (beta)
trueSiteEffects = betarnd(8,2,1,params.seqLength) + 0.1;
params.fitMtrx(2,1:params.seqLength) = trueSiteEffects;
params.trueSiteEffects = trueSiteEffects;

% Transmission scaling through logistic decline in beta
params.logisticGrowth = false;
params.logisticK = 0.5;
params.logisticR = 0.5;
params.initSize = 0.005;

% Host birth/death rates
params.mu = zeros(params.subPops,1);

% Removal rates
params.nu = ones(params.subPops,1) * 0.05;

% Migration rates
params.gamma = ones(params.subPops,params.subPops) * 0.0;

% Birth rates
randBirthRate = 0.1 + (0.2-0.1)*rand(1,1); % Uniform[0.1 - 0.2]
params.beta_0 = randBirthRate; % birth rate (lambda_0)
params.beta = diag(ones(params.subPops,1) * params.beta_0); % only on diagonal

params_file = strcat(file_out, '_params');
save(params_file, 'params');

end


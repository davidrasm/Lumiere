function [ output_args ] = main_fitBD_sim(params_file,file_out)
% Simulate trees and write tree and seqs to file

if ~exist('params_file','var')
  params_file = 'fitBD_neutralFitModel_5sites_params';
end

if ~exist('file_out','var')
  file_out = 'fitBD_neutralFitModel_5sites';
end

% Read in epi params -- moved below
% load(params_file)

% Get a tree by forward Gillespie simulation
simFlag = 1;
while(simFlag)
    
    % Moved this inside loop so won't get stuck at bad param combos
    write_randomParams_file(file_out);
    
    % Read in epi params
    load(params_file)
        
    [times, X, Y, iList, iSeqs, samples, align, simFlag] = run_stochTreeSim_seqEvo(params); %Run stochastic sim
    sampleCount = samples.sampleCount;
    sampleTimes = samples.sampleTimes;
    samplesAtTimes = samples.samplesAtTimes;
    
    if (sampleCount < params.minTreeSize || sampleCount > params.maxTreeSize)
        strout = strcat('Tree contained', {' '}, num2str(sampleCount), {' '}, 'samples');
        display(strout)
        simFlag = 1;
    end
    
end

% Plot prevalence
hmap = colormap(parula(3));
if (params.displayTrees)
    figure();
    fTime = params.finalTime / params.dt;
    for m = 1:params.subPops
        plot(times(1:fTime), X(m,1:fTime,2), '-', 'LineWidth', 2.0, 'Color', hmap(m,:)); hold on;
    end
    box off; xlabel('Time'); ylabel('Prevalence');
end

% Get entire tree with lineage states
treeParams.colorMap = hmap;

% Construct a new phylo object from iList and the samples
iList(:,5) = iList(:,5) + 1; % shift mutLoad count
phy = phylo(iList,iSeqs,sampleTimes,samplesAtTimes);
fileName = 'simTree.tre';
[phy, simTreeString] = writeSimMapTree(phy, fileName);

% Plot tree with states at mutation times
params.maxLoad = params.seqLength + 1; %max(iList(:,5));
mut_hmap = colormap(parula(params.maxLoad));
simTreeParams.colorMap = mut_hmap;
simTreeParams.colorBar = true;
if (params.displayTrees)
    simmaptree_plot(simTreeString, simTreeParams);
end

% Save mut freq to params
% params.siteMutFreqs = sum(align.Genotypes) / sampleCount;
% save(params_file, 'params');

tree_file = strcat(file_out, '_phy');
save(tree_file, 'phy');

% Save full alignment for ease of retrival
align_file = strcat(file_out, '_align');
save(align_file, 'align');

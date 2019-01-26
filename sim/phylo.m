classdef phylo
    % A class for phylogenetic trees
    % Works with MTBD and MFBD likelihood calculations
    
    properties
        nodeCount          % Number of tree nodes
        internalCount      % Number of internal nodes
        tipCount           % Number of external nodes
        nodes              % Array of tree nodes
        params             % Defines parameters for tree
    end
    
    methods
        
        % Constructor method from iList
        function obj=phylo(iList,iSeqs,sampleTimes,samplesAtTimes)
            
            % Added this so method works with or without seqs
            if isempty(iSeqs)
                pushSeqs = false;
            else
                pushSeqs = true;
            end
            
            % Initialize array for tree nodes 
            samples = 0;
            for sn = 1:length(samplesAtTimes)
                samples = samples + length(samplesAtTimes{sn});
            end
            obj.tipCount = samples;
            obj.nodes = cell(samples,1); % holds node objects
            
            % Copy iList entries for easy reference
            birthTimes = iList(:,1);
            parents = iList(:,3);
            plotState = 5; % mut load
            states = iList(:,plotState);
            
            % Build tree
            indexesInSample = cell(0); %holds active indexes in iList -- COULD THIS BE A VECTOR SO DON'T HAVE TO CONVERT?
            nodeHeights = cell(0); %hold heights of nodes
            activeTreeNodes = cell(0); %holds active indexes in tree
            activeLines = 0;
            sampleNumber = 0;
            tipTime = max(sampleTimes); %tipTime is the time of the last sample
            for n = 1:length(sampleTimes)
                
                %Find start and end times for this sampling interval
                intStartTime = sampleTimes(n);
                if (n < length(sampleTimes))
                    intEndTime = sampleTimes(n+1);
                else
                    intEndTime = 0.0;
                end
                
                %Get samples at intStartTime
                currSamples = samplesAtTimes{n};
                sampleSize = length(currSamples);
                for m = 1:sampleSize

                    indexOfSample = currSamples(m);
                    if (ismember(indexOfSample, cell2mat(indexesInSample)))
                       display('New sampled infection is already in sample!!!'); % This shouldn't happen but just to be safe
                       continue;
                    end
                    activeLines = activeLines + 1;
                    indexesInSample{activeLines} = indexOfSample;
                    sampleNumber =  sampleNumber + 1;
                    activeTreeNodes{activeLines} = sampleNumber;

                    %Add sampled lineage to node Arrays        
                    samplePop = states(indexOfSample); %iLogger.states.get(indexOfSample) + 1;
                    sampleString = strcat(num2str(samplePop), '-', num2str(indexOfSample), '_', num2str(intStartTime));
                    children = cell(0);
                    height = tipTime - intStartTime;
                    obj.nodes{sampleNumber} = treeNode(sampleString,height,intStartTime,children,samplePop);
                    
                    % Do we still need these? Or do we just need a map from
                    % activeLines to obj.nodes?
                    nodeHeights{activeLines} = tipTime - intStartTime; % still useful

                    %Add entry in lineageStates and lineageMoveTimes for new sample
                    if (pushSeqs)
                        sampleSeq = iSeqs(indexOfSample,:);
                    else
                        sampleSeq = 0;
                    end
                    obj.nodes{sampleNumber} = pushLineEvent(obj.nodes{sampleNumber}, samplePop, sampleSeq, intStartTime);
                    %lineageStates{activeLines}(1) = samplePop;
                    %lineageMoveTimes{activeLines}(1) = intStartTime;

                end
                
                nextEventTime = Inf;
                while (nextEventTime > intEndTime)

                    %Determine which event happens next among the sampled lineages
                    [nextEventTime, nextEventLineage] = getNextEventTimeInSample(indexesInSample, birthTimes);

                    if (nextEventTime > intEndTime)

                        parentOfNextEventLineage = parents(nextEventLineage); %(nextEventLineages.get(k));
                        if (parentOfNextEventLineage == -1)
                            break; %back at the first infected individuals!
                        end
                    
                        indexesInSampleMatrix = cell2mat(indexesInSample);
                        if ismember(parentOfNextEventLineage, indexesInSampleMatrix) % then this is an observed coalescent event
                            
                            parentHeight = tipTime - nextEventTime;

                            childNodeIndex = find(indexesInSampleMatrix == nextEventLineage); % better to call childActiveIndex
                            childBranchLength = parentHeight - nodeHeights{childNodeIndex}; 
                            childTreeIndex = activeTreeNodes{childNodeIndex};
                            
                            % 'Parent' is really a child of the parent node
                            parentNodeIndex = find(indexesInSampleMatrix == parentOfNextEventLineage); 
                            parentBranchLength = parentHeight - nodeHeights{parentNodeIndex};
                            parentTreeIndex = activeTreeNodes{parentNodeIndex};
                            
                            % Add parent/parent distance to childrens' node
                            newParentNode = length(obj.nodes) + 1; % add
                            obj.nodes{childTreeIndex} = setParent(obj.nodes{childTreeIndex}, newParentNode, childBranchLength);
                            obj.nodes{parentTreeIndex} = setParent(obj.nodes{parentTreeIndex}, newParentNode, parentBranchLength);
                            
                            % Add parent node
                            name = '';
                            date = tipTime - parentHeight;
                            children = [childTreeIndex, parentTreeIndex];
                            state = 0;
                            obj.nodes{newParentNode} = treeNode(name,height,date,children,state);

                            nodeHeights{parentNodeIndex} = parentHeight;
                            nodeHeights(childNodeIndex) = [];
                            activeLines = activeLines - 1;
                            activeTreeNodes{parentNodeIndex} = newParentNode; % was missing this before
                            activeTreeNodes(childNodeIndex) = []; 
                            indexesInSample(childNodeIndex) = [];

                            %Update lineageStates for parent lineage 
                            if (pushSeqs)
                                parentSeq = iSeqs(parentOfNextEventLineage,:);
                            else
                                parentSeq = 0;
                            end
                            obj.nodes{newParentNode} = pushLineEvent(obj.nodes{newParentNode}, states(parentOfNextEventLineage), parentSeq, nextEventTime);

                            % Remove child from lineage arrays
                            %lineageStates(childNodeIndex) = [];
                            %lineageMoveTimes(childNodeIndex) = [];

                        else

                            childNodeIndex = find(indexesInSampleMatrix == nextEventLineage);
                            currLineageState = states(nextEventLineage);
                            nextLineageState = states(parentOfNextEventLineage);
                            if (currLineageState ~= nextLineageState)
                                
                                treeIndex = activeTreeNodes{childNodeIndex}; % Need to know node number
                                if (pushSeqs)
                                    nextSeq = iSeqs(parentOfNextEventLineage,:);
                                else
                                    nextSeq = 0;
                                end
                                obj.nodes{treeIndex} = pushLineEvent(obj.nodes{treeIndex}, nextLineageState, nextSeq, nextEventTime);
                                %lineageStates{childNodeIndex}(end+1) = nextLineageState;
                                %lineageMoveTimes{childNodeIndex}(end+1) = nextEventTime;
                                
                            end
                            indexesInSample{childNodeIndex} = parentOfNextEventLineage;

                        end
                        
                    end
                        
                end
                
            end
            
            % Set these
            obj.nodeCount = length(obj.nodes); %number of nodes
            obj.internalCount = obj.nodeCount - obj.tipCount;
            
            % Helper functions
            function [nextEventTime, nextEventLineage] = getNextEventTimeInSample(linesInSample, birthTimes)

                lns = cell2mat(linesInSample);
                nextEventsInSample = birthTimes(lns);
                nextEventLoc = find(nextEventsInSample == max(nextEventsInSample));
                if(nextEventLoc > 1)
                    nextEventLoc = nextEventLoc(1);
                end
                nextEventTime = nextEventsInSample(nextEventLoc);
                nextEventLineage = lns(nextEventLoc); % return lineages index in iList, not index in linesInSample

            end
            
        end
            
        function obj=phyloTreeStruct(treeStruct,params) % Constructor method for Matlab's treeStruct
            
            obj.nodeCount = length(treeStruct.names); %number of nodes
            obj.internalCount = length(treeStruct.tree);
            obj.tipCount = (obj.nodeCount - obj.internalCount);
            obj.nodes = cell(0); % holds node objects
            
            % Don't need nodesVisited here???
            nodesVisited = zeros(1,obj.nodeCount);
            nodesVisited(1:obj.tipCount) = 1;
            
            nodeHeights = zeros(1,obj.nodeCount);
            
            for i = 1:obj.tipCount
                
                name = treeStruct.names{i};
                state = 0;
                
                if isfield(params,'terminalTipDate')
                    splits = regexp(name, '_', 'split');
                    if (length(splits) < 3)
                        display('No height found')
                    else
                        sampleTime = str2num(splits{end});
                        height = params.terminalTipDate - sampleTime;
                    end
                else
                    height = 0;
                end
                
                nodeHeights(i) = height;
                
                children = cell(0);
                
                obj.nodes{end+1} = treeNode(name,height,sampleTime,children,state);
                
            end
            
            for i = 1:obj.internalCount
                
                children = treeStruct.tree(i,:);
                
                % This is insured by Matlab's tree structure
                if (nodesVisited(children(1)) == 0)
                    display('Have not visited node')
                end
                
                child1 = 0;
                if isempty(children(1))
                   display('Missing a child!!!')
                else
                    child1 = children(1);
                end
                
                child2 = 0;
                if isempty(children(2))
                    display('Missing a child!!!')
                else
                    child2 = children(2);
                end
                
                childHeight = nodeHeights(child1);
                distance = treeStruct.dist(child1);
                height = childHeight + distance;
                date = params.terminalTipDate - height;
                nodeHeights(i+obj.tipCount) = height;
                nodesVisited(i+obj.tipCount) = 1;
                
                % Add parent/parent distance to childrens' node
                obj.nodes{children(1)} = setParent(obj.nodes{children(1)}, i+obj.tipCount, treeStruct.dist(child1));
                obj.nodes{children(2)} = setParent(obj.nodes{children(2)}, i+obj.tipCount, treeStruct.dist(child2));
                
                name = '';
                state = 0;
                obj.nodes{end+1} = treeNode(name,height,date,children,state);
                
            end
            %nodeTimes = params.finalTime - nodeHeights;
            
            % Walk back up tree to set heights (in distance from root)
            % If have sampling times, define root to terminal tip distance
            % Date root
            % At root date to node hieghts to get dated phylogeny
            
        end
        
        function obj=randomizeTips(obj)
            
            % Get ordered array of tip states
            tipStates = zeros(obj.tipCount,1);
            for i = 1:obj.tipCount
                tipStates(i) = obj.nodes{i}.state;
            end
            
            % Permute array of tip states
            p = randperm(obj.tipCount); % random permutation of integers 1 to n (inclusive)
            permutedStates = tipStates(p);
            
            % Assign each tip a new state
            for i = 1:obj.tipCount
                obj.nodes{i}.state = permutedStates(i);
            end
            
        end
          
        function [obj,tipNames]=getTipNames(obj)
            
            tipNames = cell(1,obj.tipCount);
            for i = 1:obj.tipCount
                name = obj.nodes{i}.name;
                tipNames{i} = name;
            end
            
        end
        
        function [obj]=trackLineage(obj,tipNode)
            
            % Track a lineage through tree from a given tipNode
            obj.nodes{tipNode} = markAsSaved(obj.nodes{tipNode});
            parentNode = obj.nodes{tipNode}.parent;
            while (~isempty(parentNode))
                obj.nodes{parentNode} = markAsSaved(obj.nodes{parentNode});
                parentNode = obj.nodes{parentNode}.parent;
            end            
            
        end
        
        function [obj]=trackLineStates(obj)
            
            nextSeq = 0; % not used here
            for i = 1:obj.nodeCount
                obj.nodes{i}.lineTimes = [];
                obj.nodes{i}.lineStates = [];
                if (obj.nodes{i}.save)
                    obj.nodes{i} = pushLineEvent(obj.nodes{i}, 2, nextSeq, obj.nodes{i}.date);
                else
                    obj.nodes{i} = pushLineEvent(obj.nodes{i}, 1, nextSeq, obj.nodes{i}.date);
                end
            end
            
        end
        
        function [obj]=convertSeqsToStates(obj,genoMap)
            
            for i = 1:obj.nodeCount
                
                totalTimes = length(obj.nodes{i}.lineTimes);         
                for t = 1:totalTimes
                    seqVec = obj.nodes{i}.lineSeqs(t,:);
                    seq = strcat(num2str(seqVec(1)),num2str(seqVec(2)));
                    [genoMap, newState] = getStateOfSeq(genoMap, seq);
                    obj.nodes{i}.lineStates(t) = newState;
                end
                
            end 
            
        end
        
        function [obj,newPhy]=prune(obj,saveList)
            
            newTipCount = length(saveList);
            for i = 1:newTipCount
                currNode = saveList(i);
                obj.nodes{currNode} = markAsSaved(obj.nodes{currNode});
            end
            
            for n = 1:obj.internalCount
                
                currNode = obj.tipCount + n;
                children = obj.nodes{currNode}.children;
                savedCount = 0;
                savedChildren = [];
                for c = 1:length(children)
                    if (obj.nodes{children(c)}.save == true)
                        savedCount = savedCount + 1;
                        savedChildren(end+1) = c;
                    end
                end
                
                if (savedCount == 0)
                    %disgard
                elseif (savedCount == 1)
                    
                    %%% THIS NEEDS TO BE FIXED FOR ROOT!!!
                    if (currNode == obj.nodeCount)
                        continue;
                    end
                    
                    singleChild = children(savedChildren);
                    parent = obj.nodes{currNode}.parent;
                    
                    % Update child's parent and parent distance
                    obj.nodes{singleChild}.parent = parent;
                    newDist = obj.nodes{singleChild}.parentDistance + obj.nodes{currNode}.parentDistance;
                    obj.nodes{singleChild}.parentDistance = newDist; 
                    
                    % Update parent's children
                    parentsChildren = obj.nodes{parent}.children;
                    rIndex = find(parentsChildren == currNode);
                    obj.nodes{parent}.children(rIndex) = singleChild;
                    
                else
                    obj.nodes{currNode} = markAsSaved(obj.nodes{currNode});
                end
                
            end
            
            % Rebuild node list and new phy tree
            newPhy = obj;
            newTipCount = 0;
            newInternalCount = 0;
            newNodes = cell(0);
            
            oldNodeLocs = [];
            
            for i = 1:obj.tipCount
                if (obj.nodes{i}.save == true)
                    newNodes{end+1} = obj.nodes{i};
                    newTipCount = newTipCount + 1;
                    oldNodeLocs(end+1) = i;
                end
            end
            
            for n = 1:obj.internalCount
                currNode = obj.tipCount + n;
                if (obj.nodes{currNode}.save == true)
                    
                    newInternalCount = newInternalCount + 1;
                    nuNode = obj.nodes{currNode};

                    children = obj.nodes{currNode}.children;
                    for c = 1:length(children)
                        
                        % Update node's children with newNodeLocs
                        nuIndex = find(oldNodeLocs == children(c));
                        nuNode.children(c) = nuIndex;
                        
                        % Update children with new parent
                        newNodes{nuIndex}.parent = newInternalCount + newTipCount;
                        
                    end
                    
                    oldNodeLocs(end+1) = currNode;
                    newNodes{end+1} = nuNode; %obj.nodes{currNode};
                    
                end
            end
           
            newNodeCount = newTipCount + newInternalCount;
            if (obj.nodeCount ~= newNodeCount)
                % Reset root
                newNodes{end}.parent = 0;
                newNodes{end}.parentDistance = 0;
            end
            
            newPhy.nodeCount = newNodeCount;
            newPhy.internalCount = newInternalCount;
            newPhy.tipCount = newTipCount;
            newPhy.nodes = newNodes;
            
        end
        
        function obj=addSequenceData(obj, align, headers)
           
            % Attach sequence data in align to each tip node
            
            for i = 1:obj.tipCount
                
                name = obj.nodes{i}.name;
                alignIndex = find(ismember(headers, name)==1);
                if isempty(alignIndex)
                    display(strcat('Could not find entry for: ', name));
                end
                seq = align(alignIndex).Sequence;
                obj.nodes{i} = addSequence(obj.nodes{i}, seq);
                
            end
            
        end
        
        function obj=mapSeqToState(obj, map)
            
            % Map sequence to numerical states
            
            for i = 1:obj.tipCount
                
                seq = obj.nodes{i}.sequence;
                L = length(seq);
                states = zeros(1,L);
                for l = 1:L
                    chr = seq(l);
                    if (isKey(map,chr))
                        states(l) = map(chr);
                    else
                        display('WARNING: Found unrecognized characters in sequence string!!!')
                        states(l) = 0;
                    end
                end
                obj.nodes{i} = addStateSeq(obj.nodes{i}, states);
                
            end
            
        end
        
        function [obj,logL]=jointReconstruction(obj, Q)
            
            % To do:
            
            maxStates = length(Q);
            
            % Get full matrix
            fullSites = length(obj.nodes{1}.stateSeq);
            fullMatrix = zeros(obj.tipCount, fullSites);
            for n = 1:obj.tipCount
                fullMatrix(n,:) = obj.nodes{n}.stateSeq;
            end
            
            % Remove sites with missing data (is this really necessary???)
            cntr = 1;
            siteMap = zeros(1,fullSites);
            for s = 1:fullSites
                pattern = fullMatrix(:,s);
                if (any(ismember(0, pattern)))
                    
                    % If treating as missing data
                    %siteMap(s) = 0; 
                    
                    % If interpolating based on consensus
                    consensus = mode(pattern(pattern ~= 0));
                    pattern(pattern == 0) = consensus;
                    siteMatrix(:,cntr) = pattern;
                    siteMap(s) = cntr; 
                    cntr = cntr+1;
                    
                else
                    siteMatrix(:,cntr) = pattern;
                    siteMap(s) = cntr; 
                    cntr = cntr+1;
                end
            end
            
            % Find unique site patterns
            [rMatrix,~,rLocs] = unique(siteMatrix','rows'); 
            rMatrix = rMatrix'; % siteMatrix = rMatrix(:,rLocs)
            
            % Get site pattern stats
            [~,sites] = size(siteMatrix);
            [~,uniqueSites] = size(rMatrix);
            nonunique = sites - uniqueSites;
            strout = strcat('Found: ', num2str(nonunique), ' non-unique site patterns');
            display(strout);
            
            C = zeros(maxStates, uniqueSites, obj.nodeCount); % optimal (conditional) reconstructed character
            L = zeros(maxStates, uniqueSites, obj.nodeCount); % likelihood of best reconstruction
            
            ancestralTypes = zeros(obj.nodeCount,uniqueSites);

            for n = 1:obj.tipCount
   
                % Set optimal character to sampled state
                sampleState = repmat(rMatrix(n,:),maxStates,1);
                C(:,:,n) = sampleState;
   
                % Set likelihoods of L(i) for all parent states i
                tx = obj.nodes{n}.parentDistance;
                
                P_ij = get_fullTransProbs(tx);
                
                % Do we need this for loop? -- can we do this with repelem()?
                for s = 1:uniqueSites
                    L(:,s,n) = P_ij(:,rMatrix(n,s));
                end
   
            end
            
            tips = obj.tipCount;
            for n = 1:obj.internalCount-1
    
                child1 = obj.nodes{n+tips}.children(1);
                child2 = obj.nodes{n+tips}.children(2);

                % Compute P_ij for all i and j
                tx = obj.nodes{n+tips}.parentDistance;
                P_ij = get_fullTransProbs(tx);

                % For each i, find j that maximizes Lx * Ly * P_ij
                for s = 1:uniqueSites
                
                    L_j_child1 = repmat(L(:,s,child1)',maxStates,1);
                    L_j_child2 = repmat(L(:,s,child2)',maxStates,1);
                    L_ij = P_ij .* L_j_child1 .* L_j_child2;
                    [tempL,tempC] = max(L_ij,[],2);

                    L(:,s,n+tips) = tempL;
                    C(:,s,n+tips) = tempC;
                
                end

            end
            
            % For root
            child1 = obj.nodes{obj.nodeCount}.children(1);
            child2 = obj.nodes{obj.nodeCount}.children(2);
            P_r = ones(maxStates, uniqueSites); % these should be the equilibrium probs
            L_r = P_r .* L(:,:,child1) .* L(:,:,child2);
            [bestL,rootType] = max(L_r);
            ancestralTypes(obj.nodeCount,:) = rootType;
            logL(1) = sum(log(bestL)); % log likelihood of ML reconstruction
            logL(2) = sum(log(bestL(bestL > 0)));
            
            infSites = sum(isinf(log(bestL)));
            if (infSites > 0)
                infLocs = find(isinf(log(bestL)));
                strout = strcat('WARNING: Log likelihood is neg inf at: ', num2str(infSites), ' sites');
                display(strout);
            end

            % Reconstruct other nodes
            for n = (obj.nodeCount-1):-1:1

                parentType = ancestralTypes(obj.nodes{n}.parent,:);
                for s = 1:uniqueSites
                    ancestralTypes(n,s) = C(parentType(s),s,n);
                end
    
            end
            
            % Reconstitute full matrix from unique sites
            if (infSites > 0)
                ancestralTypes(:,infLocs) = 0; % treat infinite sites as missing data
            end
            fullAncestralTypes = ancestralTypes(:,rLocs);
            
            % Use siteMap to get full alignment with missing data
            alignedAncestralTypes = zeros(obj.nodeCount, fullSites);
            for s = 1:fullSites
                if (siteMap(s) == 0)
                    alignedAncestralTypes(:,s) = zeros(obj.nodeCount, 1);
                else
                    alignedAncestralTypes(:,s) = fullAncestralTypes(:,siteMap(s));
                end
            end
            
            for n = 1:obj.nodeCount
                obj.nodes{n}.stateSeq = alignedAncestralTypes(n,:);
            end
            
            
                %%%% Helper subfunctions %%%%
                function [tMatrix] = get_fullTransProbs(tx)
        
                    % Returned matrix gives prob of going from i -> j in forward time
                    expQ = expm(Q * tx); % once we have procomputed expQ, do for all sites
                    tMatrix = expQ;

                end
            
        end
        
        function [obj, logL] = marginalReconstruction(obj, params, Q)
           
            % Compute the conditional likelihood of each ancestral state
            % using the Felsenstein pruning algorithm
            
            ancestralTypes = zeros(obj.nodeCount,1);

            % Data structures for reconstruction
            maxStates = params.maxLoad; %params.maxStates;
            states = 1:maxStates;
            
            % Conditional likelihoods
            L = zeros(obj.nodeCount, maxStates); % likelihood of best reconstruction
            
            for n = 1:obj.tipCount
   
                sampleState = obj.nodes{n}.state;
                L(n,sampleState) = 1.0;
   
            end
            
            tips = obj.tipCount;
            for n = 1:obj.internalCount-1
    
                child1 = obj.nodes{n+tips}.children(1);
                child2 = obj.nodes{n+tips}.children(2);
                
                tx1 = obj.nodes{child1}.parentDistance;
                tx2 = obj.nodes{child2}.parentDistance;
                
                % Compute P_ij for all i and j in fowards time
                P_ij_tx1 = get_fullTransProbs(tx1);
                P_ij_tx2 = get_fullTransProbs(tx2);
                
                L_x = P_ij_tx1 * L(child1,:)';
                L_y = P_ij_tx2 * L(child2,:)';
                L(n+tips,:) = (L_x .* L_y)';

            end
            
            % For root
            child1 = obj.nodes{obj.nodeCount}.children(1);
            child2 = obj.nodes{obj.nodeCount}.children(2);
            pi_r = ones(1,maxStates); % these should be the equilibrium probs
            L_r = pi_r .* L(child1,:) .* L(child2,:);
            [bestL,rootType] = max(L_r);
            ancestralTypes(obj.nodeCount) = rootType;
            logL = log(sum(L_r)); % log likelihood of ML reconstruction

            % Reconstruct other nodes
            for n = (obj.nodeCount-1):-1:1

                parentType = ancestralTypes(obj.nodes{n}.parent);
                
                % State probs
                tx = obj.nodes{n}.parentDistance;
                P_ij = get_fullTransProbs(tx); % could store these above
                Pr_i = L(n,:) .*  P_ij(parentType,:);
                
                % For probablistic reconstruction
                %Pr_i = Pr_i ./ sum(Pr_i);
                %sample = mnrnd(1,Pr_i)
                %type = states(sample > 0);
                
                [bestL,type] = max(Pr_i);
                ancestralTypes(n) = type;
    
            end
            
            for n = 1:obj.nodeCount
                obj.nodes{n}.lineTimes = obj.nodes{n}.date;
                obj.nodes{n}.lineStates = ancestralTypes(n);
            end
            
                %%%% Helper subfunctions %%%%
                function [tMatrix] = get_fullTransProbs(tx)
        
                    % Returned matrix gives prob of going from i -> j in forward time
                    expQ = expm(Q * tx); % once we have procomputed expQ, do for all sites
                    tMatrix = expQ;

                end
             
        end
        
        function [obj, log_Pr] = marginalGridReconstruction(obj, params, probReconst)
           
            % Compute the conditional likelihood of each ancestral state
            % using the Felsenstein pruning algorithm
            
            % If probReconst == false finds best ML reconstruction
            % If probReconst = true samples a reconstuction based on probs
            
            maxStates = params.maxLoad; %params.maxStates;
            states = 1:maxStates;
            
            % Set up Q and backQ matrices
            % Construct Q and backQ matrixes
            upVec = ones(maxStates-1, 1) * params.sigma_up;
            downVec = ones(maxStates-1, 1) * params.sigma_down;
            Q1 = diag(upVec, 1);
            Q2 = diag(downVec, -1);
            Q3 = Q1 + Q2;
            d = sum(Q3,2);
            Q = Q3 - diag(d,0); % NOTE: Q is indexed from 0 to maxStates

            Q1 = diag(upVec, -1);
            Q2 = diag(downVec, 1);
            Q3 = Q1 + Q2;
            d = sum(Q3,2);
            backQ = Q3 - diag(d,0); % NOTE: Q is indexed from 0 to maxStates
            
            ancestralTypes = zeros(obj.nodeCount,1);
            
            % Data structures for reconstruction
            % Conditional likelihoods
            %L = zeros(obj.nodeCount, maxStates); % likelihood of best reconstruction
            L = cell(0,0); % cell array data structure
            
            for n = 1:obj.tipCount
   
                sampleState = obj.nodes{n}.state;
                vecL = zeros(1,maxStates);
                vecL(sampleState) = 1.0;
                L{end+1}{1} = vecL;
                
                lineIntervals = length(obj.nodes{n}.lineTimes);
                for i = 2:lineIntervals
                    
                    tx = obj.nodes{n}.lineTimes(i-1) - obj.nodes{n}.lineTimes(i);
                    L_now = L{n}{i-1};
                    
                    % old way
                    %exp_Q_tx = get_fullTransProbs(tx);
                    %L_new = L_now * exp_Q_tx; % for backwards time
                    %L{n}{i} = L_new;
                    
                    L_new = L_now * expm(backQ * tx); % * L_now';
                    L{n}{i} = L_new;
                    
                end
   
            end
            
            tips = obj.tipCount;
            for n = 1:obj.internalCount-1
    
                child1 = obj.nodes{n+tips}.children(1);
                child2 = obj.nodes{n+tips}.children(2);
                
                %tx1 = obj.nodes{child1}.parentDistance;
                %tx2 = obj.nodes{child2}.parentDistance;
                
                nodeTime = obj.nodes{n+tips}.date; 
                tx1 = obj.nodes{child1}.lineTimes(end) - nodeTime;
                tx2 = obj.nodes{child2}.lineTimes(end) - nodeTime;
                
                % Compute P_ij for all i and j in fowards time
                P_ij_tx1 = get_fullTransProbs(tx1);
                P_ij_tx2 = get_fullTransProbs(tx2);
                
                % Prune
                L_child1 = L{child1}{end};
                L_child2 = L{child2}{end};
                
                % WHY NOT COMPUTE AS BELOW? MAKE SURE THIS IS CORRECT!!!
                L_x = P_ij_tx1 * L_child1';
                L_y = P_ij_tx2 * L_child2';
                L{n+tips}{1} = (L_x .* L_y)';
                
                lineIntervals = length(obj.nodes{n+tips}.lineTimes);
                for i = 2:lineIntervals
                    
                    tx = obj.nodes{n+tips}.lineTimes(i-1) - obj.nodes{n+tips}.lineTimes(i);
                    L_now = L{n+tips}{i-1};
                    
                    % old way 
                    %exp_Q_tx = get_fullTransProbs(tx);
                    %L_new = L_now * exp_Q_tx;
                    %L{n+tips}{i} = L_new;
                    
                    L_new = L_now * expm(backQ * tx); % * L_now';
                    L{n+tips}{i} = L_new;
                    
                end

            end
            
            % For root
            child1 = obj.nodes{obj.nodeCount}.children(1);
            child2 = obj.nodes{obj.nodeCount}.children(2);
            pi_r = ones(1,maxStates); % these should be the equilibrium probs
            %pi_r(1) = 1;
            L_r = pi_r .* L{child1}{end} .* L{child2}{end};
            logL = log(sum(L_r)); % log likelihood of ML reconstruction
            
            %if (probReconst)
            %    % For probablistic reconstruction
            %    L_r = L_r ./ sum(L_r);
            %    sample = mnrnd(1,L_r);
            %    rootType = states(sample > 0);
            %else 
            %    [bestL,rootType] = max(L_r);    
            %end
            
            % For now fix root type to reduce variance in reconstruction
            log_Pr = 0;
            rootType = obj.nodes{obj.nodeCount}.lineStates(1);
            
            ancestralTypes(obj.nodeCount) = rootType;
            %obj.nodes{obj.nodeCount}.lineStates(1) = rootType;
            
            % Reconstruct other nodes
            for n = (obj.nodeCount-1):-1:1

                parentType = ancestralTypes(obj.nodes{n}.parent);
                parentTime = obj.nodes{obj.nodes{n}.parent}.date;
                
                lineIntervals = length(obj.nodes{n}.lineTimes);
                for i = lineIntervals:-1:1
                    
                    tx = obj.nodes{n}.lineTimes(i) - parentTime;
                    P_ij = get_fullTransProbs(tx); % could store these above
                    Pr_i = L{n}{i} .*  P_ij(parentType,:);
                    
                    if (probReconst)
                        % For probablistic reconstruction
                        Pr_i = Pr_i ./ sum(Pr_i);
                        sample = mnrnd(1,Pr_i);
                        type = states(sample > 0);
                        log_Pr = log_Pr + log(Pr_i(type));
                    else 
                        [bestL,type] = max(Pr_i);    
                    end
                    
                    obj.nodes{n}.lineStates(i) = type;
                    parentType = type;
                    parentTime = obj.nodes{n}.lineTimes(i);
                    
                end
 
                % State probs
                %tx = obj.nodes{n}.parentDistance;
                %P_ij = get_fullTransProbs(tx); % could store these above
                %Pr_i = L(n,:) .*  P_ij(parentType,:);
                
                %[bestL,type] = max(Pr_i);
                ancestralTypes(n) = type;
    
            end
            
            %for n = 1:obj.nodeCount
            %    obj.nodes{n}.lineTimes = obj.nodes{n}.date;
            %    obj.nodes{n}.lineStates = ancestralTypes(n);
            %end
            
                %%%% Helper subfunctions %%%%
                function [tMatrix] = get_fullTransProbs(tx)
        
                    % Returned matrix gives prob of going from i -> j in forward time
                    expQ = expm(Q * tx); % once we have procomputed expQ, do for all sites
                    tMatrix = expQ;

                end
             
        end
        
        function [obj, log_trans_Pr, log_prop_Pr, log_curr_Pr] = marginalSiteGridReconstruction(obj, params, probReconst)
            
            % Compute the conditional likelihood of each ancestral state
            % using the Felsenstein pruning algorithm
            
            % Ancestral states are reconstructed at each state
            
            % If probReconst == false finds best ML reconstruction
            % If probReconst = true samples a reconstuction based on probs
            
            maxStates = 2; % now in terms of seq states / not mut load
            states = 1:maxStates;
            
            % Set up Q and backQ matrices
            % Construct Q and backQ matrixes
            upVec = ones(maxStates-1, 1) * params.sigma_up;
            downVec = ones(maxStates-1, 1) * params.sigma_down;
            Q1 = diag(upVec, 1);
            Q2 = diag(downVec, -1);
            Q3 = Q1 + Q2;
            d = sum(Q3,2);
            Q = Q3 - diag(d,0); % NOTE: Q is indexed from 0 to maxStates

            Q1 = diag(upVec, -1);
            Q2 = diag(downVec, 1);
            Q3 = Q1 + Q2;
            d = sum(Q3,2);
            backQ = Q3 - diag(d,0); % NOTE: Q is indexed from 0 to maxStates
            

            ancestralTypes = zeros(obj.nodeCount,1);
            sites = length(obj.nodes{1}.lineSeqs(1,:));
            ancestralSeqs = zeros(obj.nodeCount,sites);
            
            % Data structures for reconstruction
            % Conditional likelihoods
            %L = zeros(obj.nodeCount, maxStates); % likelihood of best reconstruction
            L = cell(0,0); % cell array data structure
            
            for n = 1:obj.tipCount
                
                % Fast way -- but only works for binary seqs
                sampleSeq = obj.nodes{n}.lineSeqs(1,:);
                invSampleSeq = abs(sampleSeq - 1); % zeros and ones reversed
                lineL = [invSampleSeq', sampleSeq'];
                L{end+1}{1} = lineL;

                lineIntervals = length(obj.nodes{n}.lineTimes);
                for i = 2:lineIntervals
                    tx = obj.nodes{n}.lineTimes(i-1) - obj.nodes{n}.lineTimes(i);
                    L_now = L{n}{i-1}; % this is sites x states matrix
                    L_new = L_now * expm(backQ * tx);
                    L{n}{i} = L_new;
                end
                    
            end
            
            tips = obj.tipCount;
            for n = 1:obj.internalCount-1
    
                child1 = obj.nodes{n+tips}.children(1);
                child2 = obj.nodes{n+tips}.children(2);
                
                % Old way w/o line intervals
                %tx1 = obj.nodes{child1}.parentDistance;
                %tx2 = obj.nodes{child2}.parentDistance;
                
                nodeTime = obj.nodes{n+tips}.date; 
                tx1 = obj.nodes{child1}.lineTimes(end) - nodeTime;
                tx2 = obj.nodes{child2}.lineTimes(end) - nodeTime;
                
                % Compute P_ij for all i and j in fowards time
                %P_ij_tx1 = get_fullTransProbs(tx1);
                %P_ij_tx2 = get_fullTransProbs(tx2);
                
                % Prune
                L_child1 = L{child1}{end};
                L_child2 = L{child2}{end};

                %L_x = P_ij_tx1 * L_child1';
                %L_y = P_ij_tx2 * L_child2';
                
                L_x = L_child1 * expm(backQ * tx1);
                L_y = L_child2 * expm(backQ * tx2);
                
                % L{n+tips}{1} = (L_x .* L_y)';
                
                L{n+tips}{1} = (L_x .* L_y);
                
                lineIntervals = length(obj.nodes{n+tips}.lineTimes);
                for i = 2:lineIntervals
                    
                    tx = obj.nodes{n+tips}.lineTimes(i-1) - obj.nodes{n+tips}.lineTimes(i);
                    L_now = L{n+tips}{i-1};
                    
                    % old way 
                    %exp_Q_tx = get_fullTransProbs(tx);
                    %L_new = L_now * exp_Q_tx;
                    %L{n+tips}{i} = L_new;
                    
                    L_new = L_now * expm(backQ * tx); % * L_now';
                    L{n+tips}{i} = L_new;
                    
                end

            end
            
            % For root
            child1 = obj.nodes{obj.nodeCount}.children(1);
            child2 = obj.nodes{obj.nodeCount}.children(2);
            pi_r = ones(sites,maxStates); % these should be the equilibrium probs
            L_r = pi_r .* L{child1}{end} .* L{child2}{end};
            logL = sum(log(sum(L_r,2))); % log likelihood of ML reconstruction
            
            %if (probReconst)
            %    % For probablistic reconstruction
            %    L_r = L_r ./ sum(L_r);
            %    sample = mnrnd(1,L_r);
            %    rootType = states(sample > 0);
            %else 
            %    [bestL,rootType] = max(L_r);    
            %end
            
            % For now fix root type to reduce variance in reconstruction
            %log_Pr = 0;
            log_trans_Pr = 0; % log transition probs for computing likelihood of seqs
            log_prop_Pr = 0; % log proposal prob
            log_curr_Pr = 0; % log proposal prob for current seq
            rootSeq = obj.nodes{obj.nodeCount}.lineSeqs(1,:);
            rootType = sum(rootSeq) + 1; %obj.nodes{obj.nodeCount}.lineStates(1);
            
            % These are unnecessary
            ancestralSeqs(obj.nodeCount,:) = rootSeq;
            ancestralTypes(obj.nodeCount) = rootType;

            %obj.nodes{obj.nodeCount}.lineStates(1) = rootType;
            
            % Reconstruct other nodes
            siteList = [1:sites]';
            for n = (obj.nodeCount-1):-1:1

                %parentType = ancestralTypes(obj.nodes{n}.parent);
                parentSeq = ancestralSeqs(obj.nodes{n}.parent,:);
                currParentSeq = obj.nodes{obj.nodes{n}.parent}.storedLineSeqs(1,:);
                parentTime = obj.nodes{obj.nodes{n}.parent}.date;
                
                lineIntervals = length(obj.nodes{n}.lineTimes);
                for i = lineIntervals:-1:1
                    
                    tx = obj.nodes{n}.lineTimes(i) - parentTime;
                    P_ij = get_fullTransProbs(tx); % could store these above
                    Pr_i = L{n}{i} .*  P_ij((parentSeq+1),:);
                    curr_Pr_i = L{n}{i} .*  P_ij((currParentSeq+1),:);
                    
                    if (probReconst)
                        
                        % For probablistic reconstruction
                        %Pr_i = Pr_i ./ sum(Pr_i); % old way
                        Pr_i = diag(1./sum(Pr_i,2)) * Pr_i; % divide each row by its sum
                        sample = mnrnd(1,Pr_i);
                        [~,type] = max(sample,[],2); %states(sample > 0);
                        seq = (type - 1)';
                        Pr_seq = Pr_i(sub2ind([sites, maxStates], siteList, type));
                        
                        % Only for computing hastings ratios
                        curr_Pr_i = diag(1./sum(curr_Pr_i,2)) * curr_Pr_i;
                        curr_seq = obj.nodes{n}.storedLineSeqs(i,:);
                        curr_type = curr_seq + 1;
                        Pr_curr_seq = curr_Pr_i(sub2ind([sites, maxStates], siteList, curr_type'));
                        
                        %log_Pr = log_Pr + sum(log(Pr_seq));
                        log_trans_Pr = log_trans_Pr + sum(log(P_ij(sub2ind([maxStates, maxStates], (parentSeq+1)', type))));
                        log_prop_Pr = log_prop_Pr + sum(log(Pr_seq)); % these are the normalized probs 
                        log_curr_Pr = log_curr_Pr + sum(log(Pr_curr_seq));
                        
                        if (isinf(log_prop_Pr))
                            keyboard
                        end
                        if (isinf(log_curr_Pr))
                            keyboard
                        end
                    else 
                        [bestL,type] = max(Pr_i,[],2); 
                        seq = (type - 1)';
                    end
                    
                    type = sum(seq) + 1;
                    obj.nodes{n}.lineStates(i) = type;
                    obj.nodes{n}.lineSeqs(i,:) = seq;
                    
                    parentSeq = seq;
                    currParentSeq = obj.nodes{n}.storedLineSeqs(i,:);
                    parentTime = obj.nodes{n}.lineTimes(i);
                    
                end
 
                % State probs
                %tx = obj.nodes{n}.parentDistance;
                %P_ij = get_fullTransProbs(tx); % could store these above
                %Pr_i = L(n,:) .*  P_ij(parentType,:);
                
                %[bestL,type] = max(Pr_i);
                ancestralSeqs(n,:) = seq;
                ancestralTypes(n) = type;
    
            end
            
            %for n = 1:obj.nodeCount
            %    obj.nodes{n}.lineTimes = obj.nodes{n}.date;
            %    obj.nodes{n}.lineStates = ancestralTypes(n);
            %end
            
                %%%% Helper subfunctions %%%%
                function [tMatrix] = get_fullTransProbs(tx)
        
                    % Returned matrix gives prob of going from i -> j in forward time
                    expQ = expm(Q * tx); % once we have procomputed expQ, do for all sites
                    tMatrix = expQ;

                end
             
        end
        
        function [obj] = storeLineSeqs(obj)
            
            for n = 1:1:obj.nodeCount
               [obj.nodes{n}] = storeLineData(obj.nodes{n});  
            end
            
        end
        
        function [obj,lineLoads]=mutLoad(obj,prefs,startSite,endSite)
            
            fullSites = length(obj.nodes{1}.stateSeq);
            maxPrefs = max(prefs,[],2);
            lineLoads = zeros(obj.nodeCount,1);
            
            for n = 1:obj.nodeCount
                load = 0;
                for s = startSite:endSite
                    state = obj.nodes{n}.stateSeq(s);
                    if (state ~= 0)
                        load = load + (maxPrefs(s) / prefs(s,state));
                    end
                end
                lineLoads(n) = load;
                obj.nodes{n}.state = load;
            end
            
        end
        
        function [obj,mutCounts]=countMuts(obj,wt)
            
            fullSites = length(obj.nodes{1}.stateSeq);
            mutCounts = zeros(obj.nodeCount,1);
            
            for n = 1:obj.nodeCount
                muts = 0;
                for s = 2:fullSites % ignore start codon
                    state = obj.nodes{n}.stateSeq(s);
                    if (state ~= 0 && state ~= wt(s))
                        muts = muts + 1;
                    end
                end
                mutCounts(n) = muts;
                obj.nodes{n}.state = muts;
            end
            
        end
        
        function obj=fitchParsimony(obj)
            
            % Computed with the Sankoff algorithm (see Felsenstein p15-18)
            
            % Get max state
            maxS = 0;
            for n = 1:obj.tipCount
                tipState = obj.nodes{n}.state;
                if (tipState > maxS)
                    maxS = tipState;
                end
            end
            maxS = maxS + 1; % to include zero state
            
            costs = zeros(obj.nodeCount,maxS);
            tempL = zeros(1,maxS);
            tempR = zeros(1,maxS);
            
            for n = 1:obj.tipCount
                costs(n,:) = Inf;
                tipState = obj.nodes{n}.state;
                costs(n,(tipState+1)) = 0;
            end
            
            % Postorder traversal
            for n = 1:obj.internalCount
                
                currNode = obj.tipCount + n;
                leftNode = obj.nodes{currNode}.children(1);
                rightNode = obj.nodes{currNode}.children(2);
                                
                % Compute least cost for each state i of node n
                for i = 1:maxS
                    
                    for j = 1:maxS
                        llc = costs(leftNode,j);
                        rlc = costs(rightNode,j);
                        if (i ~= j)
                            llc = llc + 1; % cost of going from i to j
                            rlc = rlc + 1; % cost of going from i to j
                        end
                        tempL(j) = llc;
                        tempR(j) = rlc;
                    end
                    costs(currNode,i) = min(tempL) + min(tempR);   
                    
                end
    
            end
            
            % Assign root state
            [~,rootIndex] = min(costs(obj.nodeCount,:));
            obj.nodes{obj.nodeCount}.state = rootIndex - 1;

            % Preorder traversal
            for n = (obj.nodeCount-1):-1:(obj.tipCount+1)
                
                currNode = n;
                parentNode = obj.nodes{currNode}.parent;
                parentState = obj.nodes{parentNode}.state + 1;
                [lcCost,lcState] = min(costs(currNode,:));
                
                if (parentState == lcState)
                    obj.nodes{currNode}.state = lcState - 1;
                else
                    parentCost = costs(currNode,parentState);
                    if (parentCost < (lcCost+1))
                        obj.nodes{currNode}.state = parentState - 1;
                    else
                        obj.nodes{currNode}.state = lcState - 1;
                    end
                end
                   
            end
            
        end
        
        function obj=addTipAnnotation(obj,label,sampleNames,sampleValues)
           
            %NOTE: this is for adding annotation from a metadate file
            
            for i = 1:obj.tipCount
                
                nameWithDate = obj.nodes{i}.name;
                splits = regexp(nameWithDate, '_', 'split');
                name = splits{1};
                metaIndex = find(ismember(sampleNames, name)==1);
                if isempty(metaIndex)
                    display(strcat('Could not find entry for: ', name));
                end
                
                val = sampleValues{metaIndex};
                obj.nodes{i} = addAnnotation(obj.nodes{i}, label, val);
                
            end
            
        end
        
        function obj=addNodeAnnotation(obj,map,label)
           
            %NOTE: this is for adding annotation based on node state
            
            for i = 1:obj.nodeCount
                
                state = obj.nodes{i}.state;
                [map,val] = getNameForState(map,state);
                
                obj.nodes{i} = addAnnotation(obj.nodes{i}, label, val);
                
            end
            
        end
        
        function obj=addNodeFitness(obj,lineLoads,discreteStates)
           
            %NOTE: this is for adding annotation based on node state
            
            label = 'fit';
            
            % Normalize lineLoads
            %lineLoads = lineLoads - min(lineLoads);
            %lineLoads = lineLoads / max(lineLoads);
            %lineLoads = ceil(lineLoads * discreteStates);
            
            edges = [1.0:0.05:1.5,2.0];
            [~,lineLoads] = histc(lineLoads,edges);
            
            for i = 1:obj.nodeCount
                
                val = num2str(lineLoads(i));
                
                obj.nodes{i} = addAnnotation(obj.nodes{i}, label, val);
                
            end
            
        end
        
        function obj=addNodeAnnotationFromSiteSeq(obj,label,site)
           
            %NOTE: this is for adding annotation based on a particular site
            % in a the siteSeq
            
            for i = 1:obj.nodeCount
                
                state = obj.nodes{i}.stateSeq(site);
                val = num2str(state);
                obj.nodes{i} = addAnnotation(obj.nodes{i}, label, val);
                
            end
            
        end
        
        function [obj]=getCoalEventStateData(obj)
            
            localEventTimes = [];
            zaEventTimes = [];
            externalEventTimes = [];
            
            for i = 1:obj.internalCount
                
                nodeIndex = i + obj.tipCount;
                date = obj.nodes{nodeIndex}.date;
                children = obj.nodes{nodeIndex}.children;
                
                child1State = obj.nodes{children(1)}.state;
                child2State = obj.nodes{children(2)}.state;
                
                if (child1State == 1 && child2State == 1) % both are local
                    
                    localEventTimes(end+1) = date;
                    
                elseif (child1State == 1) % child 1 is local   
                    
                    if (child2State == 2)
                        zaEventTimes(end+1) = date;
                    else
                        externalEventTimes(end+1) = date;
                    end
                    
                elseif (child2State == 1) % child 2 is local
                    
                    if (child1State == 2)
                        zaEventTimes(end+1) = date;
                    else
                        externalEventTimes(end+1) = date;
                    end
                    
                end
                
            end
            
            yrCenters = 1985.5:1:2013.5;
            localCounts = hist(localEventTimes,yrCenters);
            zaCounts = hist(zaEventTimes,yrCenters);
            externalCounts = hist(externalEventTimes,yrCenters);
            plot(yrCenters,localCounts,'r','LineWidth',2.0); box off;
            hold on
            plot(yrCenters,zaCounts,'g','LineWidth',2.0);
            plot(yrCenters,externalCounts,'b','LineWidth',2.0); box off;
            xlabel('Year');
            ylabel('Transmission events');
            
            
        end
        
        function [obj]=gridLineStates(obj,times)
            
            % Currently ignoring endLineTime
            
            for i = 1:obj.nodeCount
                
                nodeTime = obj.nodes{i}.date; 
                currLineTimes = obj.nodes{i}.lineTimes;
                currLineIndexes = 1:1:length(currLineTimes);
                currLineStates = obj.nodes{i}.lineStates;
                currLineSeqs = obj.nodes{i}.lineSeqs;
                
                % Lineages are missing their earliest time (time of their
                % parent node) -- should add this to constructor method???
                startLineTime = nodeTime - obj.nodes{i}.parentDistance;
                if (isempty(startLineTime))
                    startLineTime = 0.0;
                end
                endLineTime = currLineTimes(1); % latest time (also nodeTime)
                
                % Debuggin
                %i
                %startLineTime
                %endLineTime
                
                lineGridTimes = intersect(times(times <= endLineTime), times(times > startLineTime));
                if (not(ismember(endLineTime, lineGridTimes)))
                    lineGridTimes(end+1) = endLineTime;
                end
                
                % Reverse lineGridTimes so they are in reverse order
                lineGridTimes = sort(lineGridTimes,'descend');
                lineGridStates = zeros(1,length(lineGridTimes));
                lineGridSeqs = zeros(length(lineGridTimes),length(currLineSeqs(1,:)));
                
                % For each lineGridTimes, find first currLineTimes <=
                % lineGridTime
                for gt = 1:length(lineGridTimes)
                    gridTime = lineGridTimes(gt);
                    futureIndexes = currLineIndexes(currLineTimes >= gridTime);
                    stateIndex = futureIndexes(end);
                    lineGridStates(gt) = currLineStates(stateIndex);
                    lineGridSeqs(gt,:) = currLineSeqs(stateIndex,:);
                end
                
                %if (length(lineGridTimes) ~= length(lineGridStates))
                %    display('lineGridTime do not match currLineTimes')
                %end
                %if length(lineGridStates) ~= length(currLineStates)
                %    
                %end
                
                obj.nodes{i}.lineTimes = lineGridTimes;
                obj.nodes{i}.lineStates = lineGridStates;
                obj.nodes{i}.lineSeqs = lineGridSeqs;
                
            end    
            
        end
        
        function [obj,treeString]=writeSimMapTree(obj, fileName)
            
            fileID = fopen(fileName, 'w');
            treeString = '';
            
            nodeNames = cell(1,obj.tipCount);
            for i = 1:obj.tipCount
                nodeNames{i} = obj.nodes{i}.name;
            end
            
            for n = 1:obj.internalCount;
                
                currNode = obj.tipCount + n;
                daughters = obj.nodes{currNode}.children;
                nodeTime = obj.nodes{currNode}.date;
                dCount = length(daughters);
                if (dCount < 2 || dCount > 2)
                   display('WARNING: Tree is non-binary!!!') 
                end
                
                daughter1String = nodeNames{daughters(1)};
                daughter2String = nodeNames{daughters(2)};
                
                daughter1States = getLineageStateString(obj.nodes{daughters(1)}.lineStates, obj.nodes{daughters(1)}.lineTimes, nodeTime);
                daughter2States = getLineageStateString(obj.nodes{daughters(2)}.lineStates, obj.nodes{daughters(2)}.lineTimes, nodeTime);
                nodeNames{obj.tipCount+n} = strcat('(', daughter1String, ':', daughter1States, ',', daughter2String, ':', daughter2States, ')');
                
            end
            
            treeString = strcat(treeString, nodeNames{end}, ';');
            fprintf(fileID, '%s\r\n', treeString);
            fclose(fileID);
            
            % Helper functions
            function [lineString] = getLineageStateString(lineageStates, moveTimes, nodeTime)
        
                lineString = '{';
                lineageSegments = length(moveTimes);
                moveTimes(end+1) = nodeTime;
                for event = 1:lineageSegments
                    segmentLength = moveTimes(event) - moveTimes(event+1);
                    if (segmentLength == 0.0)
                        display('WARNING: Segments found with zero length');
                    end
                    segmentState = lineageStates(event);
                    lineString = strcat(lineString, num2str(segmentState), ',', num2str(segmentLength));
                    if event < lineageSegments
                        lineString = strcat(lineString, ':');
                    end
                end
                lineString = strcat(lineString, '}');
                
            end
            
        end
        
        function obj=writeAnnotatedNexus(obj, fileName)
            
            fileID = fopen(fileName, 'w');

            fprintf(fileID, '%s\r\n', '#NEXUS');
            fprintf(fileID, '%s\r\n', 'Begin trees;');
            treeString = 'tree 1 = ';
            
            nodeNames = cell(1,obj.tipCount);
            for i = 1:obj.tipCount
                nodeNames{i} = obj.nodes{i}.name;
            end
            
            for n = 1:obj.internalCount;
                
                currNode = obj.tipCount + n;
                daughters = obj.nodes{currNode}.children;
                %daughters = treeStruct.tree(n,:); % old way
                
                daughter1String = nodeNames{daughters(1)};
                daughter2String = nodeNames{daughters(2)};
                
                % Get annoation string
                daughter1Annotation = getAnnotationString(daughters(1));
                daughter2Annotation = getAnnotationString(daughters(2));
                
                % Get branch lengths
                branchLength1 = obj.nodes{daughters(1)}.parentDistance;
                branchLength2 = obj.nodes{daughters(2)}.parentDistance;
                
                nodeNames{currNode} = strcat('(', daughter1String, daughter1Annotation, ':', num2str(branchLength1), ',', daughter2String, daughter2Annotation, ':', num2str(branchLength2), ')');
            end
            treeString = strcat(treeString, nodeNames{end}, ';');
            
            fprintf(fileID, '%s\r\n', treeString);
            fprintf(fileID, '%s\r\n', 'End;');
            
            fclose(fileID);
            
            function [lineString] = getAnnotationString(line)
                
                labels = obj.nodes{line}.annotationLabels;
                vals = obj.nodes{line}.annotationValues;
                
                count = length(labels);
                lineString = '';
                if (count > 0)
                
                    % Example [&date=2011.7,year=2011]
                    lineString = '[&';
                    for antn = 1:length(labels)
                    
                        lineString = strcat(lineString, labels{antn}, '=', vals{antn});
                        if (antn < length(labels))
                            lineString = strcat(lineString, ',');
                        end
                    
                    end
                    lineString = strcat(lineString, ']');
                
                end
        
            end
            
        end
        
        function [obj]=writeNewick(obj, fileName)
            
            fileID = fopen(fileName, 'w');
            treeString = '';
            
            nodeNames = cell(1,obj.tipCount);
            for i = 1:obj.tipCount
                nodeNames{i} = obj.nodes{i}.name;
            end
            
            for n = 1:obj.internalCount;
                
                currNode = obj.tipCount + n;
                daughters = obj.nodes{currNode}.children;
                dCount = length(daughters);
                
                if (dCount < 2 || dCount > 2)
                   display('WARNING: Tree is non-binary!!!') 
                end
                
                daughter1String = nodeNames{daughters(1)};
                daughter2String = nodeNames{daughters(2)};
                
                % Get branch lengths
                branchLength1 = obj.nodes{daughters(1)}.parentDistance;
                branchLength2 = obj.nodes{daughters(2)}.parentDistance;
                
                nodeNames{obj.tipCount+n} = strcat('(', daughter1String, ':', num2str(branchLength1), ',', daughter2String, ':', num2str(branchLength2), ')');
                
            end
            
            treeString = strcat(treeString, nodeNames{end}, ';');
            fprintf(fileID, '%s\r\n', treeString);
            fclose(fileID);
            
        end
        
    end
    
end


classdef fitnessModel
    %   FitnessModel class for MFBD model
    %   Specifies fitness of each site and/or genotype
    %   Assumes 2 states at each site and multiplicative fitness effects
    
    properties
        siteEffects        % Matrix representing fitness effects of each mutation
        fitPoints          % Intervals in discretized fitness landscape
        types              % Number of possible states at each site
        sites              % Number of sites in seq
        genotypes          % Number of unique genotypes
        genotypeFitness    % Fitness values of every genotype
        genoSeqArray       % All genotypes in numerical array
    end
    
    methods
        
        function obj=fitnessModel(fitMtrx,params)
            
            % Constructor method for new instance
            
            obj.siteEffects = fitMtrx;
            [obj.types, obj.sites] = size(fitMtrx);
            obj.genotypes = obj.types^obj.sites;
            
            % General fitness model
            [obj] = enumerateGenotypes(obj); 
            [obj] = getGenotypeFitness(obj,params); % get fitness values of each genotype
            
            % Set fitness eval points on line space
            if (params.approxEProbs)
                evalPoints = 10; % set in external params?
                % obj.fitPoints = sort(linspace(0.1,1.0,evalPoints),'descend'); % gives worse approx
                minFit = min(obj.genotypeFitness);
                maxFit = max(obj.genotypeFitness);
                obj.fitPoints = sort(linspace(minFit,maxFit,evalPoints),'descend');
            else
                % If using unique genotype fitness vals
                obj.fitPoints = obj.genotypeFitness;
                %obj.fitPoints = sort(unique(obj.genotypeFitness),'descend');
            end
        
        end
        
        function [obj] = updateFitLandscape(obj,fitMtrx,params)
            
            obj.siteEffects = fitMtrx;
            
            % General fitness model
            [obj] = getGenotypeFitness(obj,params);
            
            % Set fitness eval points on line space
            if (params.approxEProbs)
                evalPoints = 10;
                % obj.fitPoints = sort(linspace(0.1,1.0,evalPoints),'descend');% gives worse approx
                minFit = min(obj.genotypeFitness);
                maxFit = max(obj.genotypeFitness);
                obj.fitPoints = sort(linspace(minFit,maxFit,evalPoints),'descend');
            else
                % Unique fitness vals
                obj.fitPoints = obj.genotypeFitness;
                %obj.fitPoints = sort(unique(obj.genotypeFitness),'descend');
            end
            
        end
        
        function [obj,transRates] = getTransitionRates(obj,params)
            
            if (params.approxEProbs)
                [obj,genoAssignments] = assignGenotypesToFitPoints(obj);
            else
                genoAssignments = 1:obj.genotypes;
            end
            
            genotypeList = 1:obj.genotypes;
            evalPoints = length(obj.fitPoints);
            transRates = zeros(evalPoints,evalPoints);
            for u = 1:evalPoints
                for v = 1:evalPoints
                    if (u ~= v)
                        
                        types_in_u = genotypeList(genoAssignments == u);
                        types_in_v = genotypeList(genoAssignments == v);
                        if (isempty(types_in_u) || isempty(types_in_v))
                            continue;
                        end
                        sumOfRates = 0;
                        for i = 1:length(types_in_u)
                            for j = 1:length(types_in_v)
                                type_i = types_in_u(i);
                                type_j = types_in_v(j);
                                genotype_i = obj.genoSeqArray(type_i,:);
                                genotype_j = obj.genoSeqArray(type_j,:);
                                absDist = sum(abs(genotype_i - genotype_j));
                                if (absDist < 2) % if genotype j can be reached by a single mutation
                                    gDist = sum(genotype_i - genotype_j);
                                    if (gDist < 0)
                                        sumOfRates = sumOfRates + params.sigma_up; % forward mutation
                                    else
                                        sumOfRates = sumOfRates + params.sigma_down; % backward mutation
                                    end
                                end
                            end
                        end
                        transRates(u,v) = sumOfRates / length(types_in_u);
                        
                    end
                end    
            end
            
        end
        
        function [obj,assignments] = assignGenotypesToFitPoints(obj)
            
            assignments = zeros(obj.genotypes,1);
            for g = 1:obj.genotypes
                [obj,assignments(g)]=getDiscreteFitPoint(obj,obj.genotypeFitness(g));
            end
            
        end
        
        function [obj,point]=getDiscreteFitPoint(obj,fitValue)
            
            % Returns discrete point in fitness space closest to input fitValue
            
            absDiffs = abs(obj.fitPoints - fitValue);
            [~,point] = min(absDiffs);
        
        end
        
        function [obj]=enumerateGenotypes(obj)
            
            % Assumes two states at each site
            
            obj.genoSeqArray = zeros(obj.genotypes,obj.sites);
            
            for i = 1:obj.sites
                
                tiles = obj.types^i;
                size = obj.genotypes/tiles;
                mutLocs = [];
                for j = 1:(tiles/2)
                    startLoc = (2*(j-1) * size) + 1;
                    endLoc = startLoc + size - 1;
                    mutLocs = [mutLocs,startLoc:endLoc];
                end
                obj.genoSeqArray(mutLocs,i) = 1;
                
            end
            
        end
        
        function [obj]=getGenotypeFitness(obj,params)
           
            % Assumes multiplicative fitness 
            %obj.genotypeFitness = zeros(obj.genotypes,1);
            %for g = 1:obj.genotypes
            %    li = sub2ind([obj.types, obj.sites], obj.genoSeqArray(g,:)+1, 1:obj.sites); % returns linear index of sites for this genotype
            %    obj.genotypeFitness(g) = prod(obj.siteEffects(li));
            %end
            
            % For epistatic fitness model
            obj.genotypeFitness(1) = params.epiFitness; % fitness of double mutant can be arbitrary 
            obj.genotypeFitness(2) = obj.siteEffects(1,1) * obj.siteEffects(2,2);
            obj.genotypeFitness(3) = obj.siteEffects(2,1) * obj.siteEffects(1,2);
            obj.genotypeFitness(4) = obj.siteEffects(1,1) * obj.siteEffects(1,2);
            
            
        end
        
        
    end
    
end


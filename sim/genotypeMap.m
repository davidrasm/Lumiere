classdef genotypeMap
    %   Map particular genotype seqs with fitness values
    %   11.07.17: Modified so double mutant has zero fitness
    
    properties
        
        siteEffects     % Matrix representing fitness effects of each mutation
        sites           % Number of sites in seq
        seqArray
        fitVals
        genotypes
        
    end
    
    methods
        
        function obj=genotypeMap(fitMtrx,params)
            
            % Constructor method for new instance
            
            obj.siteEffects = fitMtrx;
            [~,obj.sites] = size(fitMtrx);
            
            obj.genotypes = 4;

            % Get seq for each genotype
            obj.seqArray{1} = '11';
            obj.seqArray{2} = '10';
            obj.seqArray{3} = '01';
            obj.seqArray{4} = '00';
            
            % Get fitness for each genotype
            obj.fitVals(1) = fitMtrx(2,1) * fitMtrx(2,2);
            obj.fitVals(2) = fitMtrx(2,1) * fitMtrx(1,2);
            obj.fitVals(3) = fitMtrx(1,1) * fitMtrx(2,2);
            obj.fitVals(4) = fitMtrx(1,1) * fitMtrx(1,2); 
        
        end
        
        function [obj, state] = getStateOfSeq(obj, seqStr)
            
            state = find(ismember(obj.seqArray, seqStr)==1);
            
        end
        
        function [obj] = updateGenotypeMap(obj,fitMtrx,params)
            
            obj.siteEffects = fitMtrx;

            % Get fitness for each genotype
            obj.fitVals(1) = fitMtrx(2,1) * fitMtrx(2,2);
            obj.fitVals(2) = fitMtrx(2,1) * fitMtrx(1,2);
            obj.fitVals(3) = fitMtrx(1,1) * fitMtrx(2,2);
            obj.fitVals(4) = fitMtrx(1,1) * fitMtrx(1,2);
            
        end
        
        function [obj,transRates] = getTransitionRates(obj,params)
            
            % Hard-coded this for two-site model
            types = length(obj.seqArray);
            transRates = zeros(types,types);
            
            % New transRate matrix
            transRates(1,2) = params.sigma_down;
            transRates(1,3) = params.sigma_down;
            transRates(2,1) = params.sigma_up;
            transRates(2,4) = params.sigma_down;
            transRates(3,1) = params.sigma_up;
            transRates(3,4) = params.sigma_down;
            transRates(4,2) = params.sigma_up;
            transRates(4,3) = params.sigma_up;
            
        end
        
    end
    
end


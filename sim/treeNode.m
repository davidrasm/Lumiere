classdef treeNode
    % A class for nodes in a phylo tree object
    
    properties
        name             % Node name
        height           % Node age
        date             % Node date (for timed trees)
        parent           % Node parent
        parentDistance   % Distance to parent node
        children         % Array of node's children (can be more/less than two)
        state            % Node state (if defined)
        stateSeq         % If have seq of states
        annotationLabels % Cell array of annotation labels
        annotationValues % Cell array of annotation values
        save             % Indicator var used for pruning trees and tracking
        sequence         % String containing sequence data
        
        % Added these to track seq evolution through time
        lineStates
        lineSeqs
        lineTimes
        
        % Added these to store/restore line state/seq data
        storedLineSeqs;
        
    end
    
    methods
        
        function obj=treeNode(name,height,date,children,state) % Constructor for nodes
            
            obj.name = name; %number of nodes
            obj.height = height;
            obj.date = date;
            obj.children = children;
            obj.state = state;
            obj.stateSeq = 0;
            obj.annotationLabels = cell(0);
            obj.annotationValues = cell(0);
            obj.save = false; % used for pruning trees
            obj.sequence = '';
            
            % Added these to track seq evolution through time
            obj.lineStates = [];
            obj.lineSeqs = [];
            obj.lineTimes = [];
            
        end
        
        function obj=setParent(obj, parent, parentDistance)
            
            obj.parent = parent;
            obj.parentDistance = parentDistance;
            
        end
        
        function obj=pushLineEvent(obj,state,seq,time)
            
            entry = length(obj.lineStates) + 1;
            obj.lineStates(entry) = state;
            obj.lineSeqs(entry,:) = seq; 
            obj.lineTimes(entry) = time;
            
        end
        
        function obj = storeLineData(obj)
            
            obj.storedLineSeqs = obj.lineSeqs;
            
        end
        
        function obj = addSequence(obj, seq)
            
            %obj.sequence = seq; % changed this to map ebola seqs to ephylo
            obj.lineSeqs = seq;
            
        end
        
        function obj = addStateSeq(obj, seq)
            
            obj.stateSeq = seq;
            
        end
        
        function obj=addAnnotation(obj, label, val)
            
            obj.annotationLabels{end+1} = label;
            if (~isempty(val))
                obj.annotationValues{end+1} = strtrim(val);
            else
                obj.annotationValues{end+1} = 'NA';
            end
            
        end
        
        function obj=markAsSaved(obj)
            
            obj.save = true;
            
        end
        
    end
    
end


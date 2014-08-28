classdef Quail

properties(GetAccess = 'public', SetAccess = 'private')
    H;
    n;
    hist;
    pmf;
    integratedLikelihood;
    fairLikelihood;
end



methods
    function obj = Quail()
        obj.H = [];
        obj.n=0;
        obj.pmf = [];
        obj.hist = [];
        obj.integratedLikelihood =0.0;
        obj.fairLikelihood = 0.0;
    end
    
    
    
    
    function obj = createHist(obj , H)
        obj.H = H;
        obj.n = numel(H);
        obj.pmf = zeros(size(H));
        obj.hist = [];
        obj.integratedLikelihood =0.0;
        obj.fairLikelihood = 0.0;
        
        obj.hist(:,1) = unique(obj.H(:), 'rows');        % unique value
        obj.hist(:,2) = histc(obj.H(:), obj.hist(:,1));  % frequency
    end
        
    
    
    function obj = createPmf(obj)
%% Creates Probability Mass Function Using Given data
        if isempty(obj.H) ~= 1 & isempty(obj.hist) ~= 1
           
            for i=1:length(obj.hist(:,1))
                value = obj.hist(i,1);
                indexes = find(obj.H(:) == value);
                obj.pmf(indexes) = obj.hist(i,2)/obj.n;
            end
            
        end
    end
    
    
    function obj = createNormalizedPmf(obj)
        if isempty(obj.H) ~= 1 & isempty(obj.hist) ~= 1
           
            for i=1:length(obj.hist(:,1))
                value = obj.hist(i,1);
                indexes = find(obj.H(:) == value);
                obj.pmf(indexes) = obj.hist(i,2)/obj.n;
            end
            
            obj.pmf = obj.pmf./sum(obj.pmf(:));
        end
    end
    
    

    function obj = createUniformDistribution(obj, low, high, n)
%% Makes a PMF that represents a suite of hypotheses with equal p.
        i = 0:1:(n-1);
        H = (low + (high - low)) * (i/(n-1.00));
        obj = obj.createHist(H);
        obj = obj.createPmf();
    end
    
    
    function obj = getIntegratedLikelihood(obj , evidence)
        
        obj.integratedLikelihood =0.0;
        for i=1:length(obj.H)
            H = obj.H(i);
            pmfIndex = find(obj.hist(:,1) == obj.H(i));
            %disp(['H: ' num2str(H) '  pmfIndex: ' num2str(pmfIndex)]);
            
            obj = obj.getLikelihood(evidence, obj.H(i));
            obj.integratedLikelihood = obj.integratedLikelihood + (obj.fairLikelihood * obj.pmf(pmfIndex));
        end
        
    end
    
    
    function obj = getLikelihood(obj, evidence, p_of_H)
        
         obj.fairLikelihood = 0.0;
         
         count_estimated_evidence_occurs = evidence(1);
         count_estimated_evidence__not_occurs = sum(evidence(:)) - count_estimated_evidence_occurs;
         obj.fairLikelihood = (p_of_H^count_estimated_evidence_occurs) * ((1-p_of_H)^count_estimated_evidence__not_occurs);
         
    end
    
    
end
   
    
    
end

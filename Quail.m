classdef Quail

properties(GetAccess = 'public', SetAccess = 'private')

%% Actual location of the quail    
%%
    s;          % actual location
    
%% variables needed for the calculation of the location on the grid/hypothesis
%%
    H = [];    % hypothesis / grid
    r;         % rows of the grid
    c;         % columns of the grid
    % Meu = set of all estimated locations (x;y) = ( r(:) ; c(:) ) ... 
    % from where every time one location is being used in the iterations as the value of s  .............(2)

    pmf = [];
    hist = [];
    
%% Parameters for Observation Set as a Normal Distribution
%%
    N;          % Sample element number of the observation set
    sigma;      % Standard Deviation of the normal distribution
    x = [];     % observation set
    noise;      % sample noise sequence around the actual position
    K;          % covariance matrix
end



methods
    
    function obj = Quail(actual_x_location, actual_y_location)
        
        obj.s = [ actual_x_location; actual_y_location];          % set the actual location of the quail inside that grid
        
    end
    
    
    
    
    function quail = createObservationSet(quail, N, sigma)

        %% Generating Observation Set as a Normal Distribution
        %%
        quail.N = N;                                        % Sample element number of the observation set
        quail.sigma = sigma;                                % Standard Deviation of the normal distribution
        quail.K=[ quail.sigma^2 , 0 ; 0 , quail.sigma^2 ];  % covariance matrix
        quail.noise = 2*randn(2,N);                         % Create a 100-sample noise sequence with that standard deviation
        
        %Creating the observation set with the initial estimation adding some random noise; [ x = N( Meu, sigma^2 ); Normal distribution ]
        for i=1:N
            quail.x(:,i) = quail.s + quail.noise(:,i); % We are going to make the observation set by adding noise to the estimated initial state: x = s + n 
        end
        
    end
    
    
    function quail = createHypothesis(quail, rows, columns)
        %% Configuring the variables needed for the calculation of the location on the grid / Configuring variables for the hypothesis
        %%

        quail.H = ones(rows, columns);  % The initial probability of the grid ; all one
        quail.r = [ 1: 1 : rows];       % Though the grid is a 2-D space so we need to separate the Row and Column locations in two different row matrices
        quail.c = [ 1: 1 : columns];

        quail = quail.createHist();
    end
    
    
    function quail = createHist(quail)
        quail.pmf = zeros(size(quail.H));
        quail.hist = [];
        
        quail.hist(:,1) = unique(quail.H(:), 'rows');        % unique value
        quail.hist(:,2) = histc(quail.H(:), quail.hist(:,1));  % frequency
    end
    
    
    function quail = createNormalizedPmf(quail)
        if isempty(quail.H) ~= 1 & isempty(quail.hist) ~= 1
            n = numel(quail.H(:));
            
            for i=1:length(quail.hist(:,1))
                value = quail.hist(i,1);
                indexes = find(quail.H(:) == value);
                quail.pmf(indexes) = quail.hist(i,2)/n;
            end
            
            quail.pmf = quail.pmf./sum(quail.pmf(:));
        end
    end
    
    

    function quail = getLikelihood(quail, Meu)
        
    end
    
    
end
   
    
    
end

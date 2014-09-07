clear all
clc
%% Problem
%%

% A quail is hiding inside a 5x5 grid bush. We need to estimate the location of the quail.

%% Solution Equation
%%

% 's' will be the location of the quail.[ 2D location of that grid : s = [x;y] ]
% 'x' is the observation set selected to be the location of that quail
% According to Bayes's Theorem we know:

%            P( s ) * P( x|s )
% P( s| x) = --------------------                   ---------------- (1)
%                  P( x )

% where ,
% P( s ) = Prior
% P( x|s ) = likelihood
% P( s|x ) = Posterior


%% How to solve
%% 

% First we need to create random observation set
% Then we need create the initial Prior which is the pmf of that grid
% Then we need to calculate the Prior , Likelihood and Posterior
% Then the posterior becomes the new prior.

% After a specific iterations we will get lots of priors for different locations of that observation set.
% From them we can select the most likely prior to be our best estimate.


%% Using Quail Object
%%

quail = Quail(3,3);                          %% initialize with the actual location of the quail
quail = quail.createObservationSet(100, 2);  %% 100 samples using stadard deviation of 2

quail = quail.createHypothesis(5, 5);  %% initializing the grid/hypothesis
% Meu = set of all estimated locations (x;y) = ( r(:) ; c(:) ) ... 
% from where every time one location is being used in the iterations as the value of s  .............(2)

quail = quail.createNormalizedPmf(); % = P( s ) = probability mass/density function for prior distribution




%% Iterative Bayesian Process
%%

% From the Normal distribution we know that,
%                                             (x - Meu)^2
%                       1              -  --------------------
%    P(x) =  ----------------------- e        2 * sigma^2           ----------------- (3)
%              sigma * sqrt( 2 * pi)
%
% Since the mean of P( x|s ) is s, we can say that Meu = s, So we get,
%
%                                                 (x - s)^2
%                           1              -  --------------------
%    P( x|s ) =  ----------------------- e        2 * sigma^2        ----------------- (4) [ Likelihhod Calculation ]
%                  sigma * sqrt( 2 * pi)
%
% Though Meu is being determined by the iteration so we can't take all set of Meu to calculate P(x), So this is the remaining calculation we need to
% calculate otherwise:
% We know that, ∑s P( s|x ) = 1
%                                   P( s ) * P( x|s )
% So we can say from (1) that, ∑s ------------------- = 1
%                                         P(x)
%
% Which implies that , P(x) = P( s ) * P( x|s ) =  M( s| x)    ---------------------(5)


%% Implementation
%%

Pr = quail.pmf; % Initial Prior
Po = quail.pmf; % We are setting the Initial Posterior same as Initial Prior
M = 0*Pr;

for n=2:length(quail.x)
    Pr=Po;
    M=0*Pr;
    
	for i=1:length(Pr)  % For each entry in my prior table.
        for j=1:length(Pr)
            Meu=[ quail.r(i); quail.c(j)];
%% Computing Likelihood            
%% I don't understand how this calculation is similar to Eq-4.
            M(i,j) = ( 1/sqrt((2*pi)^2*det(quail.K)) ) * exp( -( quail.x(:,n)-Meu )'*inv(quail.K)*(quail.x(:,n)-Meu)/2); 
            
%% Computing Prior

            M(i,j) = M(i,j) * Pr(i,j); % Combine with prior    
        end
    end
    
	Po=M/sum(sum(M));
    [a,b]=find(Po==max(max(Po)));   % Pull out the indices at which Po achieves its max.
    sest=[ quail.r(a); quail.c(b)];             % The best estimate of the true state.
    
    disp(['Iteration: ' num2str(n) ' Estimated True State: ' num2str(quail.r(a)) ' , ' num2str( quail.c(b) )]);
end

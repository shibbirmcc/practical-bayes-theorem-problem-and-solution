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


%% Initial Data
%%

s = [1;1];          % suppose that our initial location is the first location of the grid


%% Generating Observation Set as a Normal Distribution
%%

N = 100;            % Sample element number of the observation set
sigma = 2;          % Standard Deviation of the normal distribution
n = 2*randn(2,N);   % Create a 100-sample noise sequence with that standard deviation
x = [];             % Creating the observation set with the initial estimation adding some random noise ; [ x = N( Meu, sigma^2 ); Normal distribution ]
for i=1:N
    x(:,i) = s + n(:,i); % We are going to make the observation set by adding noise to the estimated initial state: x = s + n 
end


%% Configuring the variables needed for the calculation of the location on the grid / Configuring variables for the hypothesis
%%

H = ones(5,5);          % The initial probability of the grid ; all one
r = [1:1:5];            % Though the grid is a 2-D space so we need to separate the Row and Column locations in two different row matrices
c = [1:1:5];

% Meu = set of all estimated locations (x;y) = ( r(:) ; c(:) ) ... 
% from where every time one location is being used in the iterations as the value of s  .............(2)


%% Creating the probability mass/density function for prior distribution
%%

quail = Quail();
quail = quail.createHist(H);
quail = quail.createNormalizedPmf(); % = P( s )


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
% We know that, ?s P( s|x ) = 1
%                                   P( s ) * P( x|s )
% So we can say from (1) that, ?s ------------------- = 1
%                                         P(x)
%
% Which implies that , P(x) = ?s P( s ) * P( x|s ) = ?s M( s| x)    ---------------------(5)


%% Implementation
%%

Pr = quail.pmf; % Initial Prior
Po = quail.pmf; % We are setting the Initial Posterior same as Initial Prior
M = 0*Pr;




% for i=1:length(Pr)
%     for j=1:length(Pr)
% %% Estimating true location from the observation set x which is Meu
% %% 
%         Meu = [r(i),c(j)]; % = s;
%         likelihood = 
%     end
% end


%% A solution provided by a website

K=[4,0;0,4]; % covariance matrix.

for n=2:length(x)
    Pr=Po;
    M=0*Pr;
    
	for i=1:length(Pr)  % For each entry in my prior table.
        for j=1:length(Pr)
            Meu=[r(i);c(j)];
%% Computing Likelihood            
%% I don't understand how this calculation is similar to Eq-4.
            M(i,j) = ( 1/sqrt((2*pi)^2*det(K)) ) * exp( -( x(:,n)-Meu )'*inv(K)*(x(:,n)-Meu)/2); 
            
%% Computing Prior

            M(i,j) = M(i,j) * Pr(i,j); % Combine with prior    
        end
    end
	Po=M/sum(sum(M));
    [a,b]=find(Po==max(max(Po)));   % Pull out the indices at which Po achieves its max.
    sest=[r(a);c(b)];             % The best estimate of the true state.
    
    disp(['Iteration: ' num2str(n) ' Estimated True State: ' num2str(r(a)) ' , ' num2str(c(b))]);
end





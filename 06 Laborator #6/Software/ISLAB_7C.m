function ISLAB_7C(N,sigma,lambda) 
%
% ISLAB_6A   Module that performs comparison between recursive 
%            identification procedures on an ARMAX process with 
%            constant parameters. 
%
% Inputs:	N      # simulation period (250, by default)
%               sigma  # standard deviation of PRB input 
%                        (1, by default); 
%               lambda # standard deviation of white noise 
%                        (1, by default)
%
% Outputs:      -------
%
% Explanation:	Data generated by an ARMAX process with constant 
%               parameters are employed to identify its parameters 
%               through 4 recursive identification procedures: 
%               RLS, RIV, RPEM, RPLR. 
%
% Author:   Dan Stefanoiu (*)
% Revised:  Dan Stefanoiu (*)
%           Lavinius Ioan Gliga (*)
%
% Created: April 10, 2004
% Revised: January 23, 2012
%          August 9, 2018
%
% Copyright: (*) "Politehnica" Unversity of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science
%

%
% BEGIN
% 

global FIG ;			% Figure number handler 
FIG = 1;                             

% 
% Faults preventing
% ~~~~~~~~~~~~~~~~~
if (nargin < 3)
   lambda = 1 ;
end 
if (isempty(lambda))
   lambda = 1 ;
end 
lambda = abs(lambda(1)) ; 
if (~lambda)
   lambda = 1 ; 
end 
if (nargin < 2)
   sigma = 1 ;
end 
if (isempty(sigma))
   sigma = 1 ;
end 
sigma = abs(sigma(1)) ; 
if (nargin < 1)
   N = 250 ;
end
if (isempty(N))
   N = 250 ;
end 
N = abs(fix(N(1))) ; 
if (~N)
   N = 250 ;
end 
% 
% Generating the identification data (ARMAX) 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[D,~,P] = gdata_vp(1,N,sigma,lambda,0) ; 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = version('-release'); % get version of Matlab
if (strcmp(v, '2014b') || strcmp(v, '2015a')) 
    [theta,ypred] = rarx(D,[1 1 1],'ff',1) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else % if the release of Matlab is at least 2015b
    % Estimating parameters via RLS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EstARX = recursiveARX([1 1 1]); % initialize the ARX estimator
    EstARX.ForgettingFactor = 1;   % choose the estimation method and the
                                   % factor
    theta = zeros(N, 2);
    ypred = zeros(N, 1);
    for i = 1 : N   % for each entry in the dataset
        [A, B, ypred(i)] = step(EstARX, D.y(i), D.u(i)); % identify the 
                                                         % parameters
        theta(i, 1) = A(2); % save the new parameters, to work with the
        theta(i, 2) = B(2); % rest of the program
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
figure(FIG),clf
   fig_look(FIG,1.5) ;
   subplot(311)
      plot(1:N,P.a(2)*ones(N,1),'--r', ... 
           1:N,theta(:,1),'-b') ; 
      title('Performances of Recursive Least Squares') ; 
      ylabel('AR (a)') ; 
      V = axis ; 
      axis([0 N+1 V(3:4)]) ; 
   subplot(312)
      plot(1:N,P.b(2)*ones(N,1),'--r', ... 
           1:N,theta(:,2),'-b') ; 
      ylabel('X (b)') ; 
      V = axis ; 
      axis([0 N+1 V(3:4)]) ; 
   subplot(313)
      plot(1:N,D.y-ypred,'-m') ; 
      xlabel('Normalized time') ; 
      ylabel('Pred. err.') ; 
      V = axis ; 
      axis([0 N+1 V(3:4)]) ; 
FIG=FIG+1 ; 
% 
% Estimating parameters via RIV
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[theta,ypred] = riv(D,[1 1 1]) ; 
figure(FIG),clf
   fig_look(FIG,1.5) ;
   subplot(311)
      plot(1:N,P.a(2)*ones(N,1),'--r', ... 
           1:N,theta(:,1),'-b') ; 
      title('Performances of Recursive Instrumental Variables') ; 
      ylabel('AR (a)') ; 
      V = axis ; 
      axis([0 N+1 V(3:4)]) ; 
   subplot(312)
      plot(1:N,P.b(2)*ones(N,1),'--r', ... 
           1:N,theta(:,2),'-b') ; 
      ylabel('X (b)') ; 
      V = axis ; 
      axis([0 N+1 V(3:4)]) ; 
   subplot(313)
      plot(1:N,D.y-ypred,'-m') ; 
      xlabel('Normalized time') ; 
      ylabel('Pred. err.') ; 
      V = axis ; 
      axis([0 N+1 V(3:4)]) ; 
FIG=FIG+1 ; 
% 
% Estimating parameters via RPEM
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[theta,ypred] = rpem(D,[1 1 1 0 0 1],'ff',1) ; 
figure(FIG),clf
   fig_look(FIG,1.5) ;
   subplot(411)
      plot(1:N,P.a(2)*ones(N,1),'--r', ... 
           1:N,theta(:,1),'-b') ; 
      title(['Performances of Recursive ' ... 
             'Prediction Error Minimization']) ; 
      ylabel('AR (a)') ; 
      V = axis ; 
      axis([0 N+1 V(3:4)]) ; 
   subplot(412)
      plot(1:N,P.b(2)*ones(N,1),'--r', ... 
           1:N,theta(:,2),'-b') ; 
      ylabel('X (b)') ; 
      V = axis ; 
      axis([0 N+1 V(3:4)]) ; 
   subplot(413)
      plot(1:N,P.c(2)*ones(N,1),'--r', ... 
           1:N,theta(:,3),'-b') ; 
      ylabel('MA (c)') ; 
      V = axis ; 
      axis([0 N+1 V(3:4)]) ; 
   subplot(414)
      plot(1:N,D.y-ypred,'-m') ; 
      xlabel('Normalized time') ; 
      ylabel('Pred. err.') ; 
      V = axis ; 
      axis([0 N+1 V(3:4)]) ; 
FIG=FIG+1 ; 
%
% END
%
function [a,b,lambda,magi,phii,mag,phi,f] = ISLAB_4A(at,bt,K,N,nr) 
%
% ISLAB_3a   Module that estimates the parameters of an 
%            ARX[1,1] model, with the help of 
%            Least Squares (LS) Method. 
%
% Inputs:	at  # true coefficients of AR part 
%                     ([-0.4 -0.32], by default)
%               bt  # true coefficents of X part 
%                     ([0.5 0.03], by default)
%               K   # number of frequency nodes (50, by default)
%               N   # simulation period (100, by default)
%               nr  # number of realizations (100, by default)
%
% Outputs:	a      # LS estimates of AR coefficients 
%                        (nr-by-2 matrix)
%               b      # LS estimates of X coefficents 
%                        (nr-by-2 matrix)
%               lambda # LS estimates of noise standard deviation 
%                        (nr-length vector)
%               magi   # magnitude of ideal frequency response 
%                        (noise free) (K-length vector)
%               phii   # phase of ideal frequency response 
%                        (noise free) (K-length vector)
%               mag    # magnitude of all nr estimations of 
%                        frequency response (K-by-nr matrix) 
%               phi    # phase of all nr estimations of 
%                        frequency response (K-by-nr matrix) 
%               f      # frequency axis (semi-logarithmic) 
%                        between 10^(-2) and pi (K-length vector) 
%
% Explanation:	The parameters of an ARX[2,2] model are 
%               estimated by means of LS Method. This simulator 
%               shows the estimation errors produced on 
%               Bode diagram of model, as well as the estimation 
%               of noise variance. 
%
% Author:   Dan Stefanoiu (*)
% Revised:  Dan Stefanoiu (*)
%           Lavinius Ioan Gliga (**)
%
% Last upgrade: March 8, 2004 / January 21, 2012
%               April 3, 2018
%
% Copyright: (*) "Politehnica" Unversity of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science
%

%
% BEGIN
% 

global FIG ;			% Figure number handler 
FIG = 1;    

% Constants
% ~~~~~~~~~
lambdat = 1 ;			% True variance of white noise. 
p = -0.8 ;			% Pole of input filter (AR[1]). 
Ts = 1 ; 			% Sampling period. 
% 
% Faults preventing
% ~~~~~~~~~~~~~~~~~
if (nargin < 5)
   nr = 100 ;
end
if (isempty(nr))
   nr = 100 ;
end 
nr = abs(fix(nr(1))) ; 
if (~nr)
   nr = 100 ;
end 
if (nargin < 4)
   N = 100 ;
end
if (isempty(N))
   N = 100 ;
end  
N = abs(fix(N(1))) ; 
if (~N)
   N = 100 ;
end 
if (nargin < 3)
   K = 50 ;
end 
if (isempty(K))
   K = 50 ;
end 
K = abs(fix(K(1))) ; 
if (~K)
   K = 50 ;
end  
if (nargin < 2)
   bt = [0.5 0.03] ;
end 
if (isempty(bt))
   bt = [0.5 0.03] ;
end  
bt = [0 bt(1:2)] ;
if (nargin < 1)
   at = [-0.4 -0.32] ;
end
if (isempty(at))
   at = [-0.4 -0.32] ; ;
end
at = [1 at(1:2)] ; 
a = roots(at) ; 
if (abs(a(1))>=1)
   war_err(['<ISLAB_4D>: Model unstable. ' ...
            'Reciprocal model considered instead.']) ; 
   a(1) = 1/a(1) ; 
end 
if (abs(a(2))>=1)
   war_err(['<ISLAB_4D>: Model unstable. ' ...
            'Reciprocal model considered instead.']) ; 
   a(2) = 1/a(2) ; 
end 
at = poly(a) ; 
% 
% Generating the ideal frequency response
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f = logspace(-2,pi,K) ; 	% The frequency axis. 
[magi,phii] = dbode(bt,at,Ts,f) ; 
                                % phii measured in degrees. 
% 
% Generating the Gaussian white noise
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
e = lambdat*randn(N,nr) ;  
% 
% Generating the PRB input
% ~~~~~~~~~~~~~~~~~~~~~~~~
u = sign(randn(N,nr)) ; 

% 
% Generating the nr realizations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y = filter(bt,at,u) + filter(1,at,e) ; 
% 
% Estimating the ARX parameters 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
a = zeros(nr, 1) ; 
b = zeros(nr, 1) ;
lambda = zeros(nr, 1) ;
mag = zeros(K, nr) ; 
phi = zeros(K, nr) ; 
for p=1:nr
   [ru,K] = xcov(u(:,p),'biased') ;	    % Constructing the 
   ru = ru(K>=0) ; 			    % auto-covariance of input. 
   ru = ru(1:2) ; 
   [ry,K] = xcov(y(:,p),'biased') ;	    % Constructing the 
   ry = ry(K>=0) ; 			    % auto-covariance of output. 
   ry = ry(1:3) ; 
   [ryu,K] = xcov(y(:,p),u(:,p),'biased') ; % Constructing the 
   ryu = ryu(K>=-1) ;			    % cross-covariance 
   ryu = ryu(1:4) ; 			    % output-input. 
   Rtemp = -[ryu(2) ryu(3) ; ... 
         ryu(1) ryu(2)] ; 
   R = [toeplitz(ry(1:2)) Rtemp ; ... 
        Rtemp'          toeplitz(ru)] ; 		    % Main matrix.  
   r = [-ry(2) ; -ry(3) ; ryu(3) ; ryu(4)] ;% Free vector. 
   r = R\r ; 				    % Coefficients. 
   a(p, :) = r(1:1) ; 
   b(p, :) = r(3:3) ; 
   e = [y(:,p) ... 
         [0 ; y(1:(N-1),p)] ...
         [0 ; 0 ; y(1:(N-2),p)] ... 
        -[0 ; u(1:(N-1),p)] ... 
        -[0 ; 0 ; u(1:(N-2),p)]] * [1 ; r] ;% Noise estimation. 
   lambda(p) = norm(e)/sqrt(N-4) ;% Standard deviation.  
% 
% Estimating the frequency response 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   [R,r] = dbode([0 b(end,:)],[1 a(end,:)],Ts,f) ; 
   mag(:, p) = R ; 
   phi(:, p) = r ; 
end 
% 
% Evaluating the average of parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
am = mean(a) ; 
bm = mean(b) ; 
% 
% Evaluating the average of estimation errors
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
magstd = (magi*ones(1,nr)-mag)' ; 
phistd = (phii*ones(1,nr)-phi)' ;
magm = mean(magstd)' ; 
phim = mean(phistd)' ;
% 
% Evaluating the standard deviation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
magstd = std(magstd,1)' ; 
phistd = std(phistd,1)' ;
% 
% Plotting
% ~~~~~~~~
figure(FIG),clf ; 
   fig_look(FIG,2) ; 
   subplot(211) ; 
      semilogx(f,magm,'-b', ... 
               f,magm-magstd,':r',f,magm+magstd,':r') ; 
      a = axis ; 
      axis([a(1:3) a(4)+0.1*(a(4)-a(3))]) ; 
      title(['Estimating an ARX[1,1] model ' ... 
             'by the Least Squares Method.']) ; 
      xlabel('Normalized frequency [rad/s] (log)') ; 
      ylabel('FR magnitude') ; 
      set(FIG,'DefaultTextHorizontalAlignment','left') ;
      legend('estimation error', ...
             'standard deviation tube') ; 
   subplot(212) ;
      set(FIG,'DefaultTextHorizontalAlignment','center') ;
      semilogx(f,phim,'-b', ... 
               f,phim-phistd,':r',f,phim+phistd,':r') ; 
      xlabel('Normalized frequency [rad/s] (log)') ; 
      ylabel('FR phase [deg]') ; 
      set(FIG,'DefaultTextHorizontalAlignment','left') ;
      legend('estimation error', ... 
             'standard deviation tube') ; 
FIG = FIG+1 ;
figure(FIG),clf ; 
   fig_look(FIG,2) ; 
   set(FIG,'DefaultTextHorizontalAlignment','center') ;
   plot(1:nr,lambda.^2,'-m',1:nr,ones(nr,1),'-b') ; 
   a = axis ;
   e = a(4)-a(3) ; 
   axis([0 nr+1 a(3)-0.1*3 a(4)+0.25*e]) ; 
   title(['Estimating an ARX[1,1] model ' ... 
          'by the Least Squares Method.']) ; 
   xlabel('Realization index') ; 
   ylabel('Noise variance') ; 
   set(FIG,'DefaultTextHorizontalAlignment','left') ; 
   legend('estimated','true') ; 
   text(nr/10,a(4)+0.15*e,... 
        ['        True parameters: ' ... 
          sprintf('%9.4f',[at(2:3) bt(2:3)])]) ; 
   text(nr/10,a(4)+0.05*e,... 
        ['Estimated parameters: ' ... 
          sprintf('%9.4f',[am bm])]) ; 
FIG = FIG+1 ;
%
% END
%
function ISLAB_2B(x,y,SNR)
%
% ISLAB_2B   Module that simulates how an ARMA[2,2] model 
%            equivalent to an AR[2] model plus a second noise  
%            (see Exercise #4) varies with the SNR 
%            (signal-to-noise ratio). 
%
% Inputs:	x         # real part of poles for AR model 
%                           (0.5, by default)
%           y         # imaginary part of poles for AR model
%                           (0.5, by default) 
%           SNR       # given SNR (3, by default)
%
% Outputs:	---------
%
% Explanation:	An AR(2) process, with poles specified by 
%               real part x and imaginary part y, is observed 
%               in white measurement noise. The AR(2) process 
%	            has poles in x+iy and x-iy. The disturbed output 
%               signal can be viewed as generated by an ARMA(2,2) 
%               process. The function shows how the poles and 
%               zeros of the ARMA process are varying with the SNR. 
%
% Author:   Helena Haglund (*)
% Revised:  Bjorn Wittenmark (*)
%           Dan Stefanoiu (**)
%           Lavinius Ioan Gliga (**)
%
% Last upgrade: (*)  January   3, 1997
%               (**) March     8, 2004
%                    February 10, 2012
%                    February 26, 2018
%
% Copyright: (*)  Lund Institute of Technology, SWEDEN
%                 Department of Automatic Control
%            (**) "Politehnica" Unversity of Bucharest, ROMANIA
%                 Department of Automatic Control & Computer Science
%

%
% BEGIN
% 
% Faults preventing
% ~~~~~~~~~~~~~~~~~
if nargin<3
   SNR = 3;
end
if isempty(SNR)
   SNR = 3;
end
if nargin<2
   y = 0.5;
end
if isempty(y)
   y = 0.5;
end
if nargin<1
   x = 0.5;
end
if isempty(x)
   x = 0.5;
end
% 
% Deriving the equivalent MA polynomial
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
a1=-2*x; 
a2=x^2+y^2;
A=[1 a1 a2];
lv2=1;
%
% Determine le2 (lambda_e^2)
%
Rx=[1 a1 a2;a1 1+a2 0;a2 a1 1]\[lv2;0;0];
le2=Rx(1)/SNR;
%
% Spectral factorization
%
r0=lv2+le2*( 1+a1^2+a2^2 );
r1=le2*( a1+a1*a2 );
r2=le2*a2;
R=[r0 r1 r2];
[D,leps2]=spefac(R);
poler=roots(A);
nollst=roots(D);
% 
% Plotting the poles and zeros
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subplot(211)
   tt=0:.1:6.3;			% Plot the unit circle
   plot(sin(tt),cos(tt),'k-'); 
   hold on;
   axis('square')
   set(gca,'xtick',[-1 -0.5 0 0.5 1]);
   set(gca,'ytick',[-1 -0.5 0 0.5 1]);
				% Show the poles.
   plot(real(poler),imag(poler),'rx');
   hold on;
				% Show the zeros.
   plot(real(nollst),imag(nollst),'ro');
   axis('square')
   grid on;
   title('Poles (x) and zeros (o)')
% 
% Evaluating the spectrum
% ~~~~~~~~~~~~~~~~~~~~~~~
% AR spectrum
%
[w_AR,fi_AR]=d_spektr(A,1,lv2);         
%
% Noise spectral estimation
% 
w_NOISE=w_AR;
[slask1,slask2]=size(fi_AR);
fi_NOISE=le2*ones(slask1,slask2)/(2*pi);
% 
% ARMA spectrum
%
[w_ARMA,fi_ARMA]=d_spektr(A,D,leps2);         
% 
% Plotting spectra
% ~~~~~~~~~~~~~~~~
subplot(212)
   loglog(w_AR,fi_AR,w_NOISE,fi_NOISE,w_ARMA,fi_ARMA);
   legend('AR','White noise','ARMA')
   grid on;
   title('Spectral densities');
   xlabel('\omega');
%
% END
%
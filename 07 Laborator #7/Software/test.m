%% 8.1
N1=200;
N2=800;
N3=2000;
sigma=1;
lambda=1;

[D1,~,P1] = gdata_arx(1,N1,sigma,lambda,0) ; 
[D2,~,P2] = gdata_arx(1,N2,sigma,lambda,0) ; 
[D3,~,P3] = gdata_arx(1,N3,sigma,lambda,0) ; 
FIG=1;
 EstARX = recursiveARX([1 1 1]); % initialize the ARX estimator
    EstARX.ForgettingFactor = 1;   % choose the estimation method and the
                                   % factor
    theta = zeros(N1, 2);
    ypred = zeros(N1, 1);
    for i = 1 : N1   % for each entry in the dataset
        [A, B, ypred(i)] = step(EstARX, D1.y(i), D1.u(i)); % identify the 
                                                         % parameters
        theta(i, 1) = A(2); % save the new parameters, to work with the
        theta(i, 2) = B(2); % rest of the program
    end
    %--------N1--------
figure(FIG),clf
   fig_look(FIG,1.5) ;
   subplot(311)
      plot(1:N1,P1.a(2)*ones(N1,1),'--r', ... 
           1:N1,theta(:,1),'-b') ; 
      title('Performances of Recursive Least Squares') ; 
      ylabel('AR (a)') ; 
      V = axis ; 
      axis([0 N1+1 V(3:4)]) ; 
   subplot(312)
      plot(1:N1,P1.b(2)*ones(N1,1),'--r', ... 
           1:N1,theta(:,2),'-b') ; 
      ylabel('X (b)') ; 
      V = axis ; 
      axis([0 N1+1 V(3:4)]) ; 
   subplot(313)
      plot(1:N1,D1.y-ypred,'-m') ; 
      xlabel('Normalized time') ; 
      ylabel('Pred. err.') ; 
      V = axis ; 
      axis([0 N1+1 V(3:4)]) ; 
FIG=FIG+1 ; 

%N2
figure(FIG),clf
   fig_look(FIG,1.5) ;
   subplot(311)
      plot(1:N1,P1.a(2)*ones(N1,1),'--r', ... 
           1:N1,theta(:,1),'-b') ; 
      title('Performances of Recursive Least Squares') ; 
      ylabel('AR (a)') ; 
      V = axis ; 
      axis([0 N1+1 V(3:4)]) ; 
   subplot(312)
      plot(1:N1,P1.b(2)*ones(N1,1),'--r', ... 
           1:N1,theta(:,2),'-b') ; 
      ylabel('X (b)') ; 
      V = axis ; 
      axis([0 N1+1 V(3:4)]) ; 
   subplot(313)
      plot(1:N1,D1.y-ypred,'-m') ; 
      xlabel('Normalized time') ; 
      ylabel('Pred. err.') ; 
      V = axis ; 
      axis([0 N1+1 V(3:4)]) ; 
FIG=FIG+1 ; 

%N3
figure(FIG),clf
   fig_look(FIG,1.5) ;
   subplot(311)
      plot(1:N3,P3.a(2)*ones(N3,1),'--r', ... 
           1:N3,theta(:,1),'-b') ; 
      title('Performances of Recursive Least Squares') ; 
      ylabel('AR (a)') ; 
      V = axis ; 
      axis([0 N3+1 V(3:4)]) ; 
   subplot(312)
      plot(1:N3,P3.b(2)*ones(N3,1),'--r', ... 
           1:N3,theta(:,2),'-b') ; 
      ylabel('X (b)') ; 
      V = axis ; 
      axis([0 N3+1 V(3:4)]) ; 
   subplot(313)
      plot(1:N3,D3.y-ypred,'-m') ; 
      xlabel('Normalized time') ; 
      ylabel('Pred. err.') ; 
      V = axis ; 
      axis([0 N3+1 V(3:4)]) ; 
FIG=FIG+1 ;


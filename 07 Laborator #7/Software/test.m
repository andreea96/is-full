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
    theta1 = zeros(N1, 2);
    theta2 = zeros(N2, 2);
    theta3 = zeros(N3, 2);

    ypred1 = zeros(N1, 1);
    ypred2 = zeros(N2, 1);
    ypred3 = zeros(N3, 1);
    
    for i = 1 : N1   % for each entry in the dataset
        [A, B, ypred1(i)] = step(EstARX, D1.y(i), D1.u(i)); % identify the 
                                                         % parameters
        theta1(i, 1) = A(2); % save the new parameters, to work with the
        theta1(i, 2) = B(2); % rest of the program
    end
    
     for i = 1 : N2   % for each entry in the dataset
        [A, B, ypred2(i)] = step(EstARX, D2.y(i), D2.u(i)); % identify the 
                                                         % parameters
        theta2(i, 1) = A(2); % save the new parameters, to work with the
        theta2(i, 2) = B(2); % rest of the program
     end
    
      for i = 1 : N3   % for each entry in the dataset
        [A, B, ypred3(i)] = step(EstARX, D3.y(i), D3.u(i)); % identify the 
                                                         % parameters
        theta3(i, 1) = A(2); % save the new parameters, to work with the
        theta3(i, 2) = B(2); % rest of the program
      end 
    %--------N1--------
figure(FIG),clf
   fig_look(FIG,1.5) ;
   subplot(411)
      plot(1:N1,P1.a(2)*ones(N1,1),'--r', ... 
           1:N1,theta1(:,1),'-b') ; 
      title('Performances of Recursive Least Squares') ; 
      ylabel('AR (a)') ; 
      V = axis ; 
      axis([0 N1+1 V(3:4)]) ; 
   subplot(412)
      plot(1:N1,P1.b(2)*ones(N1,1),'--r', ... 
           1:N1,theta1(:,2),'-b') ; 
      ylabel('X (b)') ; 
      V = axis ; 
      axis([0 N1+1 V(3:4)]) ; 
   subplot(413)
      plot(D1.u,'-m') ;
      ylabel('Intrari') ;
   subplot(414)
      plot(D1.y,'-m');
      ylabel('Iesiri') ;
FIG=FIG+1 ; 

%N2
figure(FIG),clf
   fig_look(FIG,1.5) ;
   subplot(411)
      plot(1:N2,P2.a(2)*ones(N2,1),'--r', ... 
           1:N2,theta2(:,1),'-b') ; 
      title('Performances of Recursive Least Squares') ; 
      ylabel('AR (a)') ; 
      V = axis ; 
      axis([0 N2+1 V(3:4)]) ; 
   subplot(412)
      plot(1:N2,P2.b(2)*ones(N2,1),'--r', ... 
           1:N2,theta2(:,2),'-b') ; 
      ylabel('X (b)') ; 
      V = axis ; 
      axis([0 N2+1 V(3:4)]) ; 
    subplot(413)
      plot(D2.u,'-m') ;
      ylabel('Intrari') ;
   subplot(414)
      plot(D2.y,'-m');
      ylabel('Iesiri') ;
FIG=FIG+1 ; 

%N3
figure(FIG),clf
   fig_look(FIG,1.5) ;
   subplot(411)
      plot(1:N3,P3.a(2)*ones(N3,1),'--r', ... 
           1:N3,theta3(:,1),'-b') ; 
      title('Performances of Recursive Least Squares') ; 
      ylabel('AR (a)') ; 
      V = axis ; 
      axis([0 N3+1 V(3:4)]) ; 
   subplot(412)
      plot(1:N3,P3.b(2)*ones(N3,1),'--r', ... 
           1:N3,theta3(:,2),'-b') ; 
      ylabel('X (b)') ;
      V = axis ; 
      axis([0 N3+1 V(3:4)]) ; 
   subplot(413)
      plot(D3.u,'-m') ;
      ylabel('Intrari') ;
   subplot(414)
      plot(D3.y,'-m');
      ylabel('Iesiri') ;
FIG=FIG+1 ;


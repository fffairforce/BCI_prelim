
function [A,C,Q,W,P_0] = create_kalman(X,Y,dt)
% X-positions + Vel
% Y-neural FR
% A-predicts pos/Vel from previous pos/Vel

A = ( X(:,2:end)*X(:,1:(end-1))' )/( X(:,1:(end-1))*X(:,1:(end-1))' );
% A(:,1) = 0;
% A(:,2) = 0;
A(:,5) = 0;
% A(1,:) = 0;
% A(2,:) = 0;
A(5,:) = 0;
A(1,1) = 1;
A(2,2) = 1;
A(1,3) = dt;
A(2,4) = dt;
A(5,5) = 1;
% C- predicts movement from neural activity
% Q-
% W- 
C = Y*X'/(X*X');
Q = (1/size(X,2))*(Y-C*X)*(Y-C*X)';
W = (1/(size(X,2)-1))*( X(:,2:end)-A*X(:,1:(end-1)) )*( X(:,2:end)-A*X(:,1:(end-1)) )';

% W(:,1) = 0;
% W(:,2) = 0;
W(:,5) = 0;
% W(1,:) = 0;
% W(2,:) = 0;
W(5,:) = 0;


P_0 = cov(X');  

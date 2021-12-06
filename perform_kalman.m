function predX = perform_kalman(Y,A,C,Q,W,P_0,X_0)

predX = zeros(size(A,1), size(Y,2));
predX(:,1) = X_0;
P = P_0;


for t = 2:size(Y,2)
    
%Calculate Kalman gain K
prior_P = A*P*A' + W;
S = C*prior_P*C' + Q;
K = prior_P*C'*inv(S);

%Use current firing rates to estimate next movement
Y_error = Y(:,t) - C*predX(:,t-1);
predX(:,t) = predX(:,t-1) + K*Y_error;

%Update movement covariance matrix
P = (eye(size(A,1)) - K*C)*prior_P;
end

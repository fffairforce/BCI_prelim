function predXX = KLpredict(XX,YY,dt)
[A,C,Q,W,P_0] = train_online_kalman(XX(1:end-1,:),YY,dt);
predXX = zeros(size(A,1), size(YY,2));
X_0 = XX(1:end-1,1);
P=P_0;
predXX(:,1) = X_0;
K = zeros(size(W,1), size(Q,1), size(YY,2));
for t = 2:size(YY,2)   
    %Calculate Kalman gain K
    prior_P = A*P*A' + W;
    S = C*prior_P*C' + Q;%state
    K(:,:,t) = prior_P*C'/S;%gain
     %Update movement covariance matrix
    P = (eye(size(A,1)) - K(:,:,t)*C)*prior_P;
end
tic
for t = 2:size(YY,2) 
    %Use current firing rates to estimate next movement
    Y_error = YY(:,t) - C*predXX(:,t-1);
    predXX(:,t) = predXX(:,t-1) + K(:,:,t)*Y_error;   
end
toc
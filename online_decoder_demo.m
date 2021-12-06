% online decoder1: 
% take in SBP

% initialize-from create kalman
% set initial states
% W
% A
% C
% Q
Pj=P_0;
%VELOCITY OUTPUT
%store previous state Xt, the position and velocity
%at time t is used to calculate the position at time t+1
%ch--channels chosen to decode
%state_dim--if task set to different states(e.g. finger
%movements)/potentially edit
ch=N;
Xtprev = Xt;
Xtoutprev = Xtout;
% Normal Update
Ct = C(1:length(ch),:);
CtQinvT = CtQinv(:,1:length(ch));
Wt = W;
At = A;
CtQinvCT = CtQinvC;

% Predict:
Xt = At*Xt;                 % prediction from previous state
prior_P = At*Pj*At'+Wt;

% Innovate/Update:
%Kt = Pj*C'/(C*Pj*C' + Q);
Kt = (eye(size(prior_P))+prior_P*CtQinvCT)\prior_P*CtQinvT;
Xt = Xt + Kt*(double(lags(chans,binLag+1)) - Ct*Xt);      % correct prediction
Pj = (eye(size(prior_P)) - Kt*Ct)*prior_P;                        % Update movement covariance matrix

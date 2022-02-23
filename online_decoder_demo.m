%% online decoder1: Base KF
% take in SBP
% SpikingBandPower: The spiking band power for all 96 channels accumulated each ms, acquired via the following procedure:
% 			1) Filtered to 300-1,000Hz using the Digital Filter Editor in the Central Software Suite (Blackrock Microsystems, LLC).
% 			2) Sampled at 2 kilo-samples per second.
% 			3) Transmitted from the Cerebus system to the xPC Target system (see methods in the main paper) at 2 kilo-samples per second.
% 			4) All samples received by the xPC within 1ms are absolute-valued then summed within each channel, contained in the SpikingBandPower field.

% testing neural signal
T = 10;
dt = 0.01;  %Time step, 10 ms
freq = 2;  %#frequency of movement
N = 384;  %# neurons 10
time_array = 0:dt:T;
mu = T/2;
sigma = 1.5;
%X - movement in 2 dims, rows 1 and 2
%Use sinusoids as movement
X = cos(2*pi*freq*(0:dt:T));
X(2,:) = sin(2*pi*freq*(0:dt:T));
% %Use gaussian distribution as 1D movement
% X=normpdf((0:dt:T),mu,sigma);
% %figure;plot((0:dt:T),y);hold on; plot((0:dt:T),4.*y)
%Use band limited noise as movement
[b,a] = butter(2, 5/(0.5/dt), 'low');
X = randn(2, length(time_array));
X = filtfilt(b,a,X')';
X = X + 0.001*randn(size(X));
X = awgn(X,10,'measured');
%rows 3 and 4, velocity of movement in 2 dims
X = [X;[diff(X,1,2),[0;0]]];
X = [X; ones(1,size(X,2))];
% X = [X;[diff(X,2),[0,0]]];
%Y - tuned firing rates of N neurons
pos_tuning = 2*pi*rand(N,1);
vel_tuning = 2*pi*rand(N,1);
tuning = [cos(pos_tuning), sin(pos_tuning), cos(vel_tuning), sin(vel_tuning)];
% tuning = [pos_tuning vel_tuning];
Y = tuning(:,1)*X(1,:) + tuning(:,2)*X(2,:) + tuning(:,3)*X(3,:) + tuning(:,4)*X(4,:);
% Y = tuning(:,1)*X(1,:) + tuning(:,2)*X(2,:); 
Y = Y + 0.1*randn(size(Y));
% use simulink for data intake?

%Inputs:        X:			Neural data of size [t, N], where t is the
%							number of binned samples and neu is the number
%							of SBP channels.
%				y:			Behavioral data of size [t, b], where t is the
%							number of binned samples and b is the
%							dimensionality of behavior. This code assumes
%							the order [position, velocity], where for
%							multidimensional filters, all positions should
%							appear prior to all velocities.
%				dt:	(optional) The size of the accumulation period,
%							in ms. This defaults to 32.
%
% initialize-from create kalman
% set initial states
% W
% A
% C
% Q
[A,C,Q,W,P_0] = train_online_kalman(X,Y,dt);
predX = zeros(size(A,1), size(Y,2));
X_0(1:5,1) = X(:,1);
predX(:,1) = X_0;
%Pj=P_0;
%VELOCITY OUTPUT
%store previous state Xt, the position and velocity
%at time t is used to calculate the position at time t+1
%state_dim--if task set to different states(e.g. finger movements)/potentially edit
% ch=N; %ch--channels chosen to decode
% Xtprev = Xt;
% Xtoutprev = Xtout;
% Normal Update
% Ct = C(1:length(ch),:);
% Wt = W;
% At = A;
P = P_0;
% % Predict:
% Xt = At*Xt;                 % prediction from previous state
% prior_P = At*Pj*At'+Wt;
% 
% % Innovate/Update:
% Kt = Pj*C'/(C*Pj*C' + Q);
% %Kt = (eye(size(prior_P))+prior_P*CtQinvCT)\prior_P*CtQinvT;
% Y_error = Y(:,t) - C*predX(:,t-1);
% Xt = Xt + Kt*Y_error;      % correct prediction   (double(lags(chans,binLag+1)) - Ct*Xt)
% Pj = (eye(size(A1,1)) - Kt*Ct)*prior_P;                        % Update movement covariance matrix
tic

for t = 2:size(Y,2)   
    %Calculate Kalman gain K
    prior_P = A*P*A' + W;
    S = C*prior_P*C' + Q;%state
    K = prior_P*C'*inv(S);%gain
    
    %Use current firing rates to estimate next movement
    Y_error = Y(:,t) - C*predX(:,t-1);
    predX(:,t) = predX(:,t-1) + K*Y_error;
    
    %Update movement covariance matrix
    P = (eye(size(A,1)) - K*C)*prior_P;
end

toc
i=1;
% %VELOCITY OUTPUT 0
% 	%current state position = position + velocity*dt (from t-1)
% 	%dt = binsize (ms)
% 	if vint_type == 0
% 		if gains(1) ~= single(0)
% 			Xt(vel_idx) = gains(1)*Xt(vel_idx)+gains(2);
% 		end
% 		Xt(pos_idx) = single(Xtprev(pos_idx)+Xt(vel_idx)*double(msCount));
% 		if gains(1) ~= single(0)
% 			Xt(intersect(pos_idx,find(Xt > gains(3))')) = gains(3);
% 			Xt(intersect(pos_idx,find(Xt < gains(4))'))  = gains(4);
% 		end
% 		Xtout = Xt;
% 	else %vint_type == 1
% 		%VELOCITY OUTPUT 1, do not reset state within KF
% 		if gains(1) ~= single(0)
% 			Xt(vel_idx) = gains(1)*Xt(vel_idx)+gains(2);
% 			Xt(intersect(pos_idx,find(Xt > gains(3))')) = gains(3);
% 			Xt(intersect(pos_idx,find(Xt < gains(4))'))  = gains(4);
% 		end
% 		Xtout(pos_idx) = single(Xtoutprev(pos_idx))+Xt(vel_idx)*double(msCount);
% 		if gains(1) ~= single(0)
% 			Xtout(intersect(pos_idx,find(Xtout > gains(3))')) = gains(3);
% 			Xtout(intersect(pos_idx,find(Xtout < gains(4))'))  = gains(4);
% 		end
% 	end
%% online decoder2:Multistate Decoder
%assume input with three class SBP
for input=1:N
    if SBP fall in Low categery
        use KF(1) get coef_H
        
    elseif SBP fall in Mid categery
        use KF(2) get coef_M

    elseif SBP fall in Hi categery
        use KF(3) get coef_L
    end

end

%% 
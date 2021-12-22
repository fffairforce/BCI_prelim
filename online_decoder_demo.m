%% online decoder1: Base KF
% take in SBP
% SpikingBandPower: The spiking band power for all 96 channels accumulated each ms, acquired via the following procedure:
% 			1) Filtered to 300-1,000Hz using the Digital Filter Editor in the Central Software Suite (Blackrock Microsystems, LLC).
% 			2) Sampled at 2 kilo-samples per second.
% 			3) Transmitted from the Cerebus system to the xPC Target system (see methods in the main paper) at 2 kilo-samples per second.
% 			4) All samples received by the xPC within 1ms are absolute-valued then summed within each channel, contained in the SpikingBandPower field.

% use simulink for data intake?

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

%state_dim--if task set to different states(e.g. finger
%movements)/potentially edit
ch=N; %ch--channels chosen to decode
Xtprev = Xt;
Xtoutprev = Xtout;
% Normal Update
Ct = C(1:length(ch),:);
Wt = W;
At = A;

% Predict:
Xt = At*Xt;                 % prediction from previous state
prior_P = At*Pj*At'+Wt;

% Innovate/Update:
Kt = Pj*C'/(C*Pj*C' + Q);
%Kt = (eye(size(prior_P))+prior_P*CtQinvCT)\prior_P*CtQinvT;
Y_error = Y(:,t) - C*predX(:,t-1);
Xt = Xt + Kt*Y_error;      % correct prediction   (double(lags(chans,binLag+1)) - Ct*Xt)
Pj = (eye(size(A1,1)) - Kt*Ct)*prior_P;                        % Update movement covariance matrix

%VELOCITY OUTPUT 0
	%current state position = position + velocity*dt (from t-1)
	%dt = binsize (ms)
	if vint_type == 0
		if gains(1) ~= single(0)
			Xt(vel_idx) = gains(1)*Xt(vel_idx)+gains(2);
		end
		Xt(pos_idx) = single(Xtprev(pos_idx)+Xt(vel_idx)*double(msCount));
		if gains(1) ~= single(0)
			Xt(intersect(pos_idx,find(Xt > gains(3))')) = gains(3);
			Xt(intersect(pos_idx,find(Xt < gains(4))'))  = gains(4);
		end
		Xtout = Xt;
	else %vint_type == 1
		%VELOCITY OUTPUT 1, do not reset state within KF
		if gains(1) ~= single(0)
			Xt(vel_idx) = gains(1)*Xt(vel_idx)+gains(2);
			Xt(intersect(pos_idx,find(Xt > gains(3))')) = gains(3);
			Xt(intersect(pos_idx,find(Xt < gains(4))'))  = gains(4);
		end
		Xtout(pos_idx) = single(Xtoutprev(pos_idx))+Xt(vel_idx)*double(msCount);
		if gains(1) ~= single(0)
			Xtout(intersect(pos_idx,find(Xtout > gains(3))')) = gains(3);
			Xtout(intersect(pos_idx,find(Xtout < gains(4))'))  = gains(4);
		end
	end
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
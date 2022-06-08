%% BCI decoder station workflow
% .rec take in(asuming SBP,sampled 100Hz)
% UDP message read and write
%data=readdata

%calculate FR

%toy data for test
T = 5; %s
dt = 0.01;  %Time step, 10 ms
freq = 2;  %#frequency of movement
N = 384;  %# neurons 10
time_array = 0:dt:T;
mu = T/2;
sigma = 1.5;
X = 5*cos(2*pi*freq*(0:dt:T));
X(2,:) = 5*sin(2*pi*freq*(0:dt:T));%change amplitude for classidfying moveemnt types
XS = 5*0.2*cos(2*pi*freq*(0:dt:T));
XS(2,:) = 5*0.2*sin(2*pi*freq*(0:dt:T));
X = [X;[diff(X,1,2),[0;0]]];
X = [X; ones(1,size(X,2));+ones(1,size(XS,2))];
XS = [XS;[diff(XS,1,2),[0;0]]];
XS = [XS; ones(1,size(XS,2));-ones(1,size(XS,2))];
pos_tuning = 2*pi*rand(N,1);
vel_tuning = 2*pi*rand(N,1);
spd_tuning = -1 + 2*rand(N,1);%uniformly distributed state tuning from [-2pi,2pi]
tuning = [cos(pos_tuning), sin(pos_tuning), cos(vel_tuning), sin(vel_tuning), spd_tuning];
Y = tuning(:,1)*X(1,:) + tuning(:,2)*X(2,:) + tuning(:,3)*X(3,:) + tuning(:,4)*X(4,:) + tuning(:,5)*X(6,:);
Y = Y + 2*randn(size(Y));
%small tuning
pos_tuning2 = 10*2*pi*rand(N,1);
vel_tuning2 = 10*2*pi*rand(N,1);
%add spped tuning(+/-) & diff amp
spd_tuning2 = -1 + 2*rand(N,1);%uniformly distributed state tuning from [-2pi,2pi]
tuning2 = [cos(pos_tuning2), sin(pos_tuning2), cos(vel_tuning2), sin(vel_tuning2), spd_tuning2];
YS = tuning2(:,1)*XS(1,:) + tuning2(:,2)*XS(2,:) + tuning2(:,3)*XS(3,:) + tuning2(:,4)*XS(4,:) + tuning2(:,5)*XS(6,:);
YS = YS + 2*randn(size(YS));
%combine X/XS in timepoint
%set movement switch onset at tp 300
XX = [X,XS];
YY = [Y,YS];
%Classifier

%training KF parameters

% decode in realtime-assuming KF training done
load("KF_para_test.mat")
sw=301;
predX = zeros(size(A,1), size(YY(:,1:sw),2));
X_0(1:5,1) = XX(1:5,1);
P=P_0;
predX(:,1) = X_0;
K = zeros(size(W,1), size(Q,1), size(YY(:,1:sw),2));
neurFreq_counter = 0;
time_array = zeros(10000,1);
next_tme_targ = dt; %s
segment = 0; % Initialize segment variable, keeps track of where in trial setup the monkey is
tic;
Keepgoing = 1;
while Keepgoing
    curr_time =toc;
    if (curr_time > next_tme_targ)
        %Count the number of time steps
        missed_bins = floor((curr_time - next_tme_targ)./dt);
        if missed_bins>-1
        next_targ = next_targ + (missed_bins+1)*time_step;  %The next time target calculated, always in whole time steps
        end
        
        %For debugging
        if missed_bins>2
            display(num2str(segment))
        end
%         rawNeuronFreqs = recieveSBP('read_.rec');
%         neuronFreqs = (rawNeuronFreqs-unitMeans);
        curr_time_idx = missed_bins;%this or 1,2,3...      
        neuronFreqs = YY(:,curr_time_idx);
        neurFreq_counter = neurFreq_counter+1;
        time_array(neurFreq_counter) = toc;
        
        %calculate score - Calculate Kalman gain K
        %movement state classifier then feed in different KF parameters
        prior_P = A*P*A' + W;
        S = C*prior_P*C' + Q;%state
        K(:,:,t) = prior_P*C'*inv(S);%gain
         %Update movement covariance matrix
        P = (eye(size(A,1)) - K(:,:,t)*C)*prior_P;
        %Use current firing rates to estimate next movement
        YY_error = YY(:,t) - C*predX(:,t-1);
        predX(:,t) = predX(:,t-1) + K(:,:,t)*YY_error;%estimated position/velocity

        %store score
        score_matrix(neurFreq_counter,:) = predX;
    end
end
% output cursor position 
predPos = predX(1:2,1);
destinationAddress = '10.52.14.13';
destinationPort = 64713;
u=udpport;
write(u,predPos,destinationAddress,destinationPort)
data = read(u,u.NumBytesWritten,"uint8")

%legacy version to test
ipA = '192.168.1.204'; portA = 3030;
ipB = '10.52.14.13';  portB = 3031;  % Modify these values to be those of your second computer.
ipC = '192.168.1.250'; portC = 3033;
%%Create UDP Object
% udpA = udp(ipB,portB,'LocalPort',portA);
udpA = udp(ipA,portA,'LocalPort',portC);
udpC = udp(ipC,portC,'LocalPort',portA);

%%Connect to UDP Object
fopen(udpA)
fopen(udpC)

fprintf(udpC,'This is test message number one.')
fprintf(udpA,'This is test message number two.')
fprintf(udpC,'doremifasolatido')
fscanf(udpC)
fscanf(udpA)

uBroadcaster = udpport("datagram")
uBroadcaster.EnableBroadcast = true;
uReceiver1 = udpport("byte","LocalPort",2020,"EnablePortSharing",true)
write(uBroadcaster,1:5,"uint8","",2020);%10.31.79.255 10.52.14.255
uReceiver1Count = uReceiver1.NumBytesAvailable
data1 = read(uReceiver1,uReceiver1Count,"uint8")
%% functions to use

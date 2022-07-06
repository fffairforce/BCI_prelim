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
X(2,:) = 5*sin(2*pi*freq*(0:dt:T));%change amplitude for classifying movement types
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
% feature selection1 - LDA
response = XX(6,:);
MdlLinear = fitcdiscr(YY',XX(6,:),'DiscrimType','linear'); 

rng(1)
cvp=cvpartition(size(YY,2),'HoldOut',0.1);
idxTrn=training(cvp);
idxTest=test(cvp);
tblTrn = array2table(YY(:,idxTrn)');
tblTrn.Y = response(idxTrn)';
tblTest = array2table(YY(:,idxTest)');
tblTest.Y = response(idxTest)';
Md1 = fitcdiscr(tblTrn,'Y');
predLabel = predict(Md1,YY(:,idxTest)');
confusionchart(response(idxTest),predLabel)
% feature selection2 - SVM
Md2 = fitcsvm(tblTrn,'Y');
predLabel2 = predict(Md2,YY(:,idxTest)');
confusionchart(response(idxTest),predLabel2)

predictNew = predict(Md1,YY(:,800:1000)');
%HMM-depennding on problem to be solve, or else linear classifier/svm would
%be enough to classify movement type
state_num = 2;
%set hight initial tans matrix so that state is likely stay steady
trans= [0.95,0.05;
         0.05,0.95];
emmis = rand(state_num,length(response));%unkown observation number for each input data(use LDA observations?)
emmis = emmis./sum(emmis,state_num);
%generate random sequence for hmm to start
[seq,states] = hmmgenerate(length(MdlLinear.Y),trans,emmis); %need to take known Markov Model parameters
[TRANS_hat, EMIS_hat, LL] = hmmtrain(seq, trans, emmis,'Algorithm','BaumWelch','maxiterations',50);
estimatesStates = hmmviterbi(seq,trans,emmis);
%training KF parameters

% decode in realtime-assuming KF training done
load("KF_para_test.mat")
sw=301;
% SBP=YY;%[n t]
predX = zeros(size(A,1), size(YY(:,1:sw),2));
X_0(1:5,1) = XX(1:5,1);
P=P_0;
predX(:,1) = X_0;
K = zeros(size(W,1), size(Q,1), size(YY(:,1:sw),2));
neurFreq_counter = 0;
time_array = zeros(10000,1);
next_time_targ = dt; %s
segment = 0; % Initialize segment variable, keeps track of where in trial setup the monkey is
tic;
Keepgoing = 1;
while Keepgoing
    curr_time =toc;
    if (curr_time > next_time_targ)
        %Count the number of time steps
        missed_bins = floor((curr_time - next_time_targ)./dt);
        if missed_bins>-1
        next_time_targ = next_time_targ + (missed_bins+1)*dt;  %The next time target calculated, always in whole time steps
        end 
        
        %For debugging
        if missed_bins>2
            display(num2str(segment))
        end
%         rawNeuronFreqs = recieveSBP('read_.rec');
%         neuronFreqs = (rawNeuronFreqs-unitMeans);%[n t]
%         curr_time_idx = missed_bins;%this or 1,2,3...      
        curr_time_idx = curr_time_idx+1;%for debug
        neuronFreqs = YY(:,curr_time_idx);%will miss indexes between SBP reading/what is the time interval for SBP
        neurFreq_counter = neurFreq_counter+1;
        time_array(neurFreq_counter) = toc;
%         tic
        %calculate score - Calculate Kalman gain K
        %movement state classifier then feed in different KF parameters
        prior_P = A*P*A' + W;
        S = C*prior_P*C' + Q;%state
        K(:,:,curr_time_idx) = prior_P*C'*inv(S);%gain
         %Update movement covariance matrix
        P = (eye(size(A,1)) - K(:,:,curr_time_idx)*C)*prior_P;
        %Use current firing rates to estimate next movement
        YY_error = neuronFreqs - C*predX(:,curr_time_idx-1);
        predX(:,curr_time_idx) = predX(:,curr_time_idx-1) + K(:,:,curr_time_idx)*YY_error;%estimated position/velocity
%         toc
        %store score
        score_matrix(neurFreq_counter,:) = predX(:,curr_time_idx);
        %send UDP
        predPos = predX(1:2,curr_time_idx);

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
ipA = '10.31.75.149'; portA = 3030;
ipB = '10.52.14.10';  portB = 3031;  % Modify these values to be those of your second computer.
ipC = '192.168.1.250'; portC = 3033;
%%Create UDP Object
% udpC = udp(ipA,portA,'LocalPort',portC);
% udpA = udp(ipC,portC,'LocalPort',portA);
% udpB = udp(ipB,portB,'LocalPort',portA)
% %%Connect to UDP Object
% fopen(udpA)
% fopen(udpC)
% fopen(udpB)
% fprintf(udpA,'This is test message number one.')
% fprintf(udpA,'This is test message number two.')
% fprintf(udpB,'doremifasolatido')
% fscanf(udpA)
% fscanf(udpC)
% fscanf(udpB)%working for lab P--> WL laptop


uBroadcaster = udpport("LocalHost",ipA,"LocalPort",portA)
uBroadcaster.EnableBroadcast = true;
write(uBroadcaster,predPos,"double",ipA,portB)

uReceiver1 = udpport("byte","LocalHost",ipA,"LocalPort",portA+1,"EnablePortSharing",true)
% write(uBroadcaster,1:5,"uint8","",2020);%10.31.79.255 10.52.14.255
uReceiver1Count = uReceiver1.NumBytesAvailable
data1 = read(uReceiver1,uReceiver1Count,"double")%B->A
data1 = reshape(data1,2,[])
%% functions to use

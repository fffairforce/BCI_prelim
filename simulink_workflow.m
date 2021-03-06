% UDP flow in simulink
% take inputs from Trodes
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

%input = YY(:,1);
%% if test KinArm, input Pos data constently 
Keepgoing = 1;
tic
%network config
% ipA = '10.31.75.185'; portA = 3030;
% ipB = '10.52.14.10';  portB = 3031;  % Modify these values to be those of your second computer.
% %actual udp sending
% uBroadcaster = udpport("LocalHost",ipA,"LocalPort",portA);
% uBroadcaster.EnableBroadcast = true;

while Keepgoing
    curr_time =toc;
    if curr_time < length(XX)
        curr_time_idx = ceil(curr_time);
    else
        curr_time_idx = ceil(mod(curr_time,length(XX)));
    end
    input = [XX(1,curr_time_idx),XX(2,curr_time_idx)];
    var = sim("SendUdp2Kinarm_new.slx");
    display(input)
    % send out UDP from PC1

    % Decoder receive data in PC2
    % Note:set samplerate in 'SendUdp2Kinarm_new/UDP Receive' to meet inpput
    % samplerate
%     outdata_PC1
    
    % Decode position prepare PredPos
    
    % send out UDP from PC2
    write(uBroadcaster,input,"double",ipA,portB)
    
    sendPos = outdata_PC1;

    % KinArm station recieve data in PC3 
    % replace handpos in 'COTPerturb/PNDL_PerturbHandInTarget/Embedded MATLAB InsideTarget' with PredPos
    % dislink handfeedback block from Kinarm
    % start with no perturbu
    %actual udp recieving
    %might add a sleep period for the udp to send over
%     Receiver1 = udpport("byte","LocalHost",ipA,"LocalPort",portB,"EnablePortSharing",true);
%     % write(uBroadcaster,1:5,"uint8","",2020);%10.31.79.255 10.52.14.255
%     Receiver1Count = Receiver1.NumBytesAvailable;
%     data1 = read(Receiver1,Receiver1Count,"double");%B->A
%     data1 = reshape(data1,2,[]);
%     print('data received!')

end
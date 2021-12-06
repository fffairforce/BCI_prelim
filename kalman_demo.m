T = 10;
dt = 0.01;  %Time step, 10 ms
freq = 2;  %#frequency of movement
N = 10;  %# neurons
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

figure;plot(Y')

%%Create Kalman
[A,C,Q,W,P_0] = create_kalman(X,Y,dt);


X_0(1:5,1) = X(:,1);
%%Perform Kalman filtering
predX = perform_kalman(Y,A,C,Q,W,P_0,X_0);


figure
plot(X(1,:))
hold on
plot(predX(1,:))
legend('Actual', 'Predicted')

figure
plot(X(2,:))
hold on
plot(predX(2,:))
legend('Actual', 'Predicted')
%% try PCA on ideal model 


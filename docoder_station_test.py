import numpy as np
import time

#calculate FR
if __name__ == "__main__":
    T = 5.00
    dt = 0.01
    freq = 2
    N = 384
    time_array = np.arange(0.00,T,dt)
    mu = T/2
    sigma = 1.5
    X = 5*np.cos(2*np.pi*freq*time_array)
    X = np.stack((X,5*np.sin(2*np.pi*freq*time_array)))
    XS = 5*0.2*np.cos(2*np.pi*freq*time_array)
    XS = np.stack((X,5*0.2*np.sin(2*np.pi*freq*time_array)))
    X = np.concatenate((np.diff(X,axis=1).T,np.zeros((2,1)).T),axis=0).T
    X = np.concatenate((X, np.ones((1,np.size(X,1))),+np.ones((1,np.size(XS,1)))))
    XS = np.concatenate((np.diff(XS,axis=1).T,np.zeros((2,1)).T),axis=0).T
    XS = np.concatenate((XS, np.ones((1,np.size(XS,1))),-np.ones((1,np.size(XS,1)))))
    pos_tuning = 2*pi*rand(N,1);
    vel_tuning = 2*pi*rand(N,1);
    spd_tuning = -1 + 2*rand(N,1);%uniformly distributed state tuning from [-2pi,2pi]
    tuning = [cos(pos_tuning), sin(pos_tuning), cos(vel_tuning), sin(vel_tuning), spd_tuning];
    Y = tuning(:,1)*X(1,:) + tuning(:,2)*X(2,:) + tuning(:,3)*X(3,:) + tuning(:,4)*X(4,:) + tuning(:,5)*X(6,:);
    Y = Y + 2*randn(size(Y));
    # small tuning
    pos_tuning2 = 10*2*pi*rand(N,1);
    vel_tuning2 = 10*2*pi*rand(N,1);
    %add spped tuning(+/-) & diff amp
    spd_tuning2 = -1 + 2*rand(N,1);%uniformly distributed state tuning from [-2pi,2pi]
    tuning2 = [cos(pos_tuning2), sin(pos_tuning2), cos(vel_tuning2), sin(vel_tuning2), spd_tuning2];
    YS = tuning2(:,1)*XS(1,:) + tuning2(:,2)*XS(2,:) + tuning2(:,3)*XS(3,:) + tuning2(:,4)*XS(4,:) + tuning2(:,5)*XS(6,:);
    YS = YS + 2*randn(size(YS));
    # combine X/XS in timepoint
    # set movement switch onset at tp 300
    XX = [X,XS];
    YY = [Y,YS];
    [A,C,P_0,Q,W] = np.load("KF_para_test.mat")

    sw=301
    # SBP=YY;%[n t]
    predX = np.zeros(np.size(A,1), np.size(YY(:,1:sw)))#size(A,1), size(YY(:,1:sw),2)
    X_0[1:5,1 ]= XX[1:5,1]
    P=P_0
    predX[:,1] = X_0
    K = np.zeros(np.size(W,1), np.size(Q,1), np.size(YY[:,1:sw],2))
    neurFreq_counter = 0;
    time_array = np.zeros(10000,1)
    next_time_targ = dt
    segment = 0 # Initialize segment variable, keeps track of where in trial setup the monkey is
    time_start = time.time()
    Keepgoing = 1

    while Keepgoing:
        curr_time = time.time() - time_start
        if (curr_time > next_time_targ):
            #Count the number of time steps
            missed_bins = np.floor((curr_time - next_time_targ)./dt)
            if missed_bins>-1:
                next_time_targ = next_time_targ + (missed_bins+1)*dt  #The next time target calculated, always in whole time steps

            #For debugging
            if missed_bins>2:
                display(np.num2str(segment))

             # rawNeuronFreqs = recieveSBP('read_.rec');
             # neuronFreqs = (rawNeuronFreqs-unitMeans);%[n t]
             # curr_time_idx = missed_bins;%this or 1,2,3...
            curr_time_idx = curr_time_idx+1 #for debug
            neuronFreqs = YY[:,curr_time_idx] #will miss indexes between SBP reading/what is the time interval for SBP
            neurFreq_counter = neurFreq_counter+1
            time_array[neurFreq_counter] = time.time()
            # calculate score - Calculate Kalman gain K
            # movement state classifier then feed in different KF parameters
            prior_P = A*P*np.transpose(A) + W;
            S = C*prior_P*np.transpose(C) + Q #state
            K[:,:,curr_time_idx] = prior_P*np.transpose(C)*np.linalog.inv(S) #gain
            # Update movement covariance matrix
            P = (np.eye(np.size(A,1)) - K[:,:,curr_time_idx]*C)*prior_P
            # Use current firing rates to estimate next movement
            YY_error = neuronFreqs - C*predX[:,curr_time_idx-1]
            predX[:,curr_time_idx] = predX[:,curr_time_idx-1] + K[:,:,curr_time_idx]*YY_error # estimated position/velocity
            # store score
            score_matrix[neurFreq_counter,:] = predX[:,curr_time_idx]
            #send UDP
            predPos = predX[1:2,curr_time_idx]


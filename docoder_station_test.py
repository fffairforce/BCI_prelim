import numpy as np
import time
from scipy.io import loadmat, savemat
import socket

#calculate FR
if __name__ == "__main__":
    T = 5.00
    dt = 0.01
    freq = 2
    N = 384
    time_array = np.arange(0.00, T + dt, dt)
    mu = T / 2
    sigma = 1.5
    X = 5 * np.cos(2 * np.pi * freq * time_array)
    X = np.stack((X, 5 * np.sin(2 * np.pi * freq * time_array)))
    XS = 5 * 0.2 * np.cos(2 * np.pi * freq * time_array)
    XS = np.stack((XS, 5 * 0.2 * np.sin(2 * np.pi * freq * time_array)))
    X = np.concatenate((X, np.concatenate((np.diff(X, axis=1).T, np.zeros((2, 1)).T), axis=0).T), axis=0)
    X = np.concatenate((X, np.ones((1, np.size(X, 1))), +np.ones((1, np.size(XS, 1)))), axis=0)
    XS = np.concatenate((XS, np.concatenate((np.diff(XS, axis=1).T, np.zeros((2, 1)).T), axis=0).T), axis=0)
    XS = np.concatenate((XS, np.ones((1, np.size(XS, 1))), -np.ones((1, np.size(XS, 1)))), axis=0)
    np.random.seed(10)
    pos_tuning = 2 * np.pi * np.random.rand(N, 1)
    vel_tuning = 2 * np.pi * np.random.rand(N, 1)
    spd_tuning = -1 + 2 * np.random.rand(N, 1) # uniformly distributed state tuning from [-2pi,2pi]
    tuning = np.concatenate(
        (np.cos(pos_tuning), np.sin(pos_tuning), np.cos(vel_tuning), np.sin(vel_tuning), spd_tuning), axis=1)
    Y = np.matmul(tuning[:, 0:1], X[0:1, :]) + np.matmul(tuning[:, 1:2], X[1:2, :]) + np.matmul(tuning[:, 2:3],X[2:3, :]) \
        + np.matmul(tuning[:, 3:4], X[3:4, :]) + np.matmul(tuning[:, 4:], X[5:, :])
    Y = Y + 2 * np.random.randn(Y.shape[0], Y.shape[1])  # is there a easier expression?
    # small tuning
    np.random.seed(20)
    pos_tuning2 = 10 * 2 * np.pi * np.random.rand(N, 1)
    vel_tuning2 = 10 * 2 * np.pi * np.random.rand(N, 1)
    # add spped tuning(+/-) & diff amp
    spd_tuning2 = -1 + 2 * np.random.rand(N, 1)  # uniformly distributed state tuning from [-2pi,2pi]
    tuning2 = np.concatenate(
        (np.cos(pos_tuning2), np.sin(pos_tuning2), np.cos(vel_tuning2), np.sin(vel_tuning2), spd_tuning2), axis=1)
    YS = np.matmul(tuning2[:, 0:1], XS[0:1, :]) + np.matmul(tuning2[:, 1:2], XS[1:2, :]) + np.matmul(tuning2[:, 2:3],XS[2:3,:])\
         + np.matmul(tuning2[:, 3:4], XS[3:4, :]) + np.matmul(tuning2[:, 4:], XS[5:, :])
    YS = YS + 2 * np.random.randn(YS.shape[0], YS.shape[1])  # is there a easier expression?
    # combine X/XS in timepoint
    # set movement switch onset at tp 300
    XX = np.concatenate((X, XS), axis=1)
    YY = np.concatenate((Y, YS), axis=1)
    # classifier for movement types

    #load KF trained parameters(for also submovements)
    KF = loadmat("KF_para_test.mat")
    A, C, P_0, Q, W = KF['A'], KF['C'], KF['P_0'], KF['Q'], KF['W']

    curr_time_idx = 2
    YY_error = []
    score_matrix = predX[:, :1]  # initilize
    sw=301
    # SBP=YY;%[n t]
    predX = np.zeros((np.size(A, 0), np.size(YY[:, :sw], 1)))  # size(A,1), size(YY(:,1:sw),2)
    X_0 = XX[:5, :1]
    P = P_0
    predX[:, :1] = X_0
    K = np.zeros((np.size(W, 1), np.size(Q, 1), np.size(YY[:, :sw], 1)))
    neurFreq_counter = 0
    time_array = np.zeros((10000,1))
    next_time_targ = dt
    segment = 0 # Initialize segment variable, keeps track of where in trial setup the monkey is
    time_start = time.time()
    Keepgoing = 1

    while Keepgoing:
        curr_time = time.time() - time_start
        if curr_time > next_time_targ:
            # Count the number of time steps
            missed_bins = np.floor((curr_time - next_time_targ)./dt)
            if missed_bins>-1:
                next_time_targ = next_time_targ + (missed_bins+1)*dt  # The next time target calculated, always in whole time steps

            # For debugging
            if missed_bins>2:
                print(np.num2str(segment))

             # rawNeuronFreqs = recieveSBP('read_.rec');
             # neuronFreqs = (rawNeuronFreqs-unitMeans);%[n t]
             # curr_time_idx = missed_bins;%this or 1,2,3...
            curr_time_idx = curr_time_idx+1 #for debug
            neuronFreqs = YY[:,
                          curr_time_idx:curr_time_idx + 1]  # will miss indexes between SBP reading/what is the time interval for SBP
            neurFreq_counter = neurFreq_counter + 1
            time_array[neurFreq_counter] = time.time()
            # calculate score - Calculate Kalman gain K
            # movement state classifier then feed in different KF parameters
            prior_P = np.matmul(np.matmul(A, P), A.T) + W
            S = np.matmul(np.matmul(C, prior_P), C.T) + Q  # state
            K[:, :, curr_time_idx - 1:curr_time_idx] = np.matmul(np.matmul(prior_P, C.T), np.linalg.inv(S)).reshape(
                np.size(K, 0), np.size(K, 1), 1)  # gain
            # Update movement covariance matrix
            P = np.eye(np.size(A, 1)) - np.matmul(np.matmul(K[:, :, curr_time_idx], C), prior_P)
            # Use current firing rates to estimate next movement
            YY_error = neuronFreqs - np.matmul(C, predX[:, curr_time_idx - 1:curr_time_idx])
            predX[:, curr_time_idx:curr_time_idx + 1] = predX[:, curr_time_idx - 1:curr_time_idx] + np.matmul(
                np.squeeze(K[:, :, curr_time_idx:curr_time_idx + 1], axis=2), YY_error)  # estimated position/velocity
            # store score
            score_matrix = np.concatenate((score_matrix, predX[:, curr_time_idx:curr_time_idx + 1]), axis=1)
            # send UDP
            predPos = predX[0:2, curr_time_idx:curr_time_idx + 1]

    # UDP receive in py
    UDP_IP = "127.0.0.1"
    UDP_PORT = 5005

    sock = socket.socket(socket.AF_INET, # Internet
                         socket.SOCK_DGRAM) # UDP
    sock.bind((UDP_IP, UDP_PORT))

    while True:
        data, addr = sock.recvfrom(2048) # buffer size is 1024 bytes
        print("received message: %s" % data)
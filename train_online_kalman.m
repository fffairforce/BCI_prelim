% create Kalman parameters
function [A,C,Q,W,P_0] = create_kalman(X,Y,binsize)
% This function trains our online SBP Kalman filter.
%	Inputs:		X:			Behavioral data of size [t, b], where t is the
%							number of binned samples and b is the
%							dimensionality of behavior. This code assumes
%							the order [position, velocity], where for
%							multidimensional filters, all positions should
%							appear prior to all velocities.
%				Y:			Neural data of size [t, n], where t is the
%							number of binned samples and n is the number
%							of SBP channels. //TBD
%				binsize:	(optional) The size of the accumulation period,
%							in ms. This defaults to 32.
%	Outputs:	Aout:		The trained trajectory model parameters. This
%							matrix is size [11, 11] populated from the top
%							left according to the size of the state vector.
%				Wout:		The trained trajectory model uncertainty
%							parameters. This matrix is size [11, 11]
%							populated from the top left according to the
%							size of the state vector.
%				Cout:		The trained observation model parameters. This
%							matrix is size [96, 11], columns populated from
%							the left according to the size of the state
%							vector and from the top according to goodChans.
%				CtQinvout:	A helpful input matrix to save online
%							operations, computed as C'/Q, where C is the
%							observation model parameters and Q is the
%							uncertainty in the observation model. This
%							matrix is size [11, 96], populated from the
%							left according to goodChans and from the top
%							according to the size of the state vector.
%				CtQinvCout:	A helpful input matrix to save online
%							operations, computed as C'/Q*C, where C is the
%							observation model parameters and Q is the
%							uncertainty in the observation model. This
%							matrix is size [11, 11], populated from the top
%							left according to the size of the state vector.
%				BinLagout:	The quantity of bins to lag the neural data for
%							optimal offline decoding, implemented online.
%				r:			The correlation coefficient between the true
%							and predicted behavior for the optimal Kalman
%							filter.
%				Q:			Observation model noise covariance.


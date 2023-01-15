#!/usr/bin/env python
# coding: utf-8

import numpy as np
from math import radians, pi
from scipy.fftpack import fft, fftfreq
from statsmodels.tsa.stattools import acf

# in 3d
def norm(vec):
    return(np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2))

def dotProduct(vec1,vec2):
    return(vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2])

def crossProduct(vec1,vec2):
    vec=[0,0,0]
    vec[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1]
    vec[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2]
    vec[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0]
    return vec

# rearrange the initial data to the same size of fiducial data
def dataCut(trange,data):
    stateVec = [0]*6
    for i in range(6):
        stateVec[i]=np.interp(trange,data[:,0],data[:,i+1])
    stateVec=np.transpose(stateVec)
#     stateVec.insert(0,initalVec)
    return(stateVec)


# estimate the period and the average value of a period
def period(values,dt,top_k_seasons=5):
    
    fft_series = fft(values)
    power = np.abs(fft_series)
    sample_freq = fftfreq(fft_series.size)

    pos_mask = np.where(sample_freq > 0)
    freqs = sample_freq[pos_mask]
    powers = power[pos_mask]
    
    top_k_idxs = np.argpartition(powers, -top_k_seasons)[-top_k_seasons:]
    top_k_power = powers[top_k_idxs]
    fft_periods = (1 / freqs[top_k_idxs]).astype(int)

#     print(f"top_k_power: {top_k_power}")
#     print(f"fft_periods: {fft_periods}")
    scores = []
    lags = []
    
    for lag in fft_periods:
#         lag = fft_periods[np.abs(fft_periods - time_lag).argmin()]
        acf_score = acf(values, nlags=lag)[-1]
        scores.append(acf_score)
        lags.append(lag)
#         print(f"lag: {lag} fft acf: {acf_score}")
    point = scores/np.sum(scores)
    return(np.abs(round(point@lags,5))*dt)


# calculate the altitude
def altitude(stateVec):
    Pos = stateVec[0:3]
    Altitude=(np.sqrt(Pos[0]**2+Pos[1]**2+Pos[2]**2))
    return(Altitude-6371000)

# calculate the difference of altitude
def errorAltitude(stateVec0,stateVec1):
    return(altitude(stateVec0)-altitude(stateVec1))


# calculate the distance between two statevectors
def distance(stateVec0,stateVec1):
    Pos0 = stateVec0[0:3]
    Pos1 = stateVec1[0:3]
    Dr = [Pos0[i]-Pos1[i] for i in range(0,len(Pos1))]
    return(norm(Dr))


# calculate the long-track difference(ATD),cross-track difference(CTD) and radial difference(RD)
# in the reference of stateVec1 (predicted orbit)
def errorEstimator(stateVec0,stateVec1):
    error = {'ATD':0,'CTD':0,'RD':0}

    Pos0 = stateVec0[0:3]
    V0 = stateVec0[3:6]
    Pos1 = stateVec1[0:3]
    V1 = stateVec1[3:6]
    
    n = crossProduct(Pos1,V1)
    Dr = [Pos0[i]-Pos1[i] for i in range(0,len(Pos1))]
    
    error['ATD'] = dotProduct(Dr,V1)/norm(V1)
    error['CTD'] = dotProduct(Dr,Pos1)/norm(Pos1)
    error['RD'] = dotProduct(Dr,n)/norm(n)
    
    return error


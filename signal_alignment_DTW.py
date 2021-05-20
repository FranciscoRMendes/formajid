# -*- coding: utf-8 -*-
"""
Created on Wed May 19 20:45:14 2021

@author: Francisco
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 08:14:58 2020
@author: frmendes@deloitte.com
"""

#DO NOT CIRCULATE

'''
Here we walk through an example of why dynamic time warping is needed in order to do the sequential analysis. 
While all the coughs are definitely by the same person, they are all occurring at naturally different speeds. 
the DTW will allow us to study the differnce and then align the coughs as necessary
'''

def smoother(sample, lf=1600, 
             var_thres=0.00003, ss_thres=0.00003,
             mode='var'):
    assert mode in ['ss', 'var']
    
    truth = sample.copy()
    len_sample = len(truth)
    l = int((len(truth))/lf + 1)
    
    
    if mode=='var':
        for i in range(l):
            if i*lf >= len_sample:
                break
                
            variance = np.var(truth[i*lf : i*lf + lf])
            if variance > var_thres:
                truth[i*lf : min(i*lf + lf, len_sample)] = 1
            else:
                truth[i*lf : min(i*lf + lf, len_sample)] = 0
    
    if mode=='ss':
        for i in range(l):
            if i*lf >= len_sample:
                break
                
            sum_squares = sum(abs(truth[i*lf:i*lf+lf]) ** 2)
            if sum_squares > ss_thres * lf:
                truth[i*lf : min(i*lf + lf, len_sample)] = 1
            else:
                truth[i*lf : min(i*lf + lf, len_sample)] = 0
    print(i*lf + lf, end=' | ')
    return truth

import pandas as pd
import scipy
import librosa
import sounddevice as sd
#import soundfile as sf
import numpy as np
import os
import matplotlib.pyplot as plt
os.chdir('C:/Users/frmendes/Documents/Merck Cough Accelerator/dat')
from os import listdir
from hmmlearn import hmm
from os.path import isfile, join
onlyfiles = [f for f in listdir() if isfile(join(f))]
from librosa import display
an = pd.read_excel('../annotation.xlsx')
fs = 44000



'''
Here, we load a "gold standard cough". This is the best cough as it is very representative of what Allison's coughs 
sound like and we use this as a starting point for any reference for alignment and phase etc
'''


i = 7
f = onlyfiles[i]
y,fs = librosa.load(f,fs)

#--------------------------------------------------------------#
y = scipy.signal.filtfilt([1,-1],[1 , -0.98],y);
y_norm = y / np.sqrt(sum(y**2)/len(y))
y_cough = y[smoother(y_norm, lf=1600, mode='var', var_thres=0.075)==1]
y_cough = y_cough[an.start.loc[i]:an.end.loc[i]]
#--------------------------------------------------------------#

librosa.display.waveplot(y_cough)
gold_y = y_cough
gold_mfcc = librosa.feature.mfcc(gold_y,sr=fs,n_mfcc=21)



'''
Load up any other cough
'''

i = 6
f = onlyfiles[i]
y,fs = librosa.load(f,fs)

#--------------------------------------------------------------#
y = scipy.signal.filtfilt([1,-1],[1 , -0.98],y);
y_norm = y / np.sqrt(sum(y**2)/len(y))
y_cough = y[smoother(y_norm, lf=1600, mode='var', var_thres=0.075)==1]
y_cough = y_cough[an.start.loc[i]:an.end.loc[i]]
#--------------------------------------------------------------#

librosa.display.waveplot(y_cough)





'''
Here we extract the chroma features, we can choose any features we like
edit : actually extracts mfcc but didnt change var names
'''

n_fft = 2048
hop_size = 1024

x_1_chroma = librosa.feature.mfcc(y=gold_y, sr=fs,n_fft=n_fft,hop_length=hop_size)
x_2_chroma = librosa.feature.mfcc(y=y_cough, sr=fs,n_fft = n_fft,hop_length=hop_size)

plt.figure(figsize=(16, 8))
plt.subplot(2, 1, 1)
plt.title('Mel Representation of $X_1$')
librosa.display.specshow(x_1_chroma, x_axis='time',
                         y_axis='mel', cmap='gray_r',hop_length=hop_size)
plt.colorbar()
plt.subplot(2, 1, 2)
plt.title('Mel Representation of $X_2$')
librosa.display.specshow(x_2_chroma, x_axis='time',
                         y_axis='mel', cmap='gray_r',hop_length=hop_size)
plt.colorbar()
plt.tight_layout()

'''
The Variable best_cost contains the distance between the two coughs, lower is better
'''

D, wp = librosa.sequence.dtw(X=x_1_chroma, Y=x_2_chroma, metric='cosine')
best_cost = D[wp[-1, 0], wp[-1, 1]]

wp_s = np.asarray(wp) * hop_size / fs

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
librosa.display.specshow(D, x_axis='time', y_axis='time',
                         cmap='gray_r', hop_length=hop_size)
imax = ax.imshow(D, cmap=plt.get_cmap('gray_r'),
                 origin='lower', interpolation='nearest', aspect='auto')
ax.plot(wp_s[:, 1], wp_s[:, 0], marker='o', color='r')
plt.title('Warping Path on Acc. Cost Matrix $D$')
plt.colorbar()





import matplotlib
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8))

# Plot x_1
librosa.display.waveplot(gold_y, sr=fs, ax=ax1)
ax1.set(title='Gold Standard Cough')

# Plot x_2
librosa.display.waveplot(y_cough, sr=fs, ax=ax2)
ax2.set(title='Sample Cough')

plt.tight_layout()

trans_figure = fig.transFigure.inverted()
lines = []
arrows = 30
points_idx = np.int16(np.round(np.linspace(0, wp.shape[0] - 1, arrows)))

# for tp1, tp2 in zip((wp[points_idx, 0]) * hop_size, (wp[points_idx, 1]) * hop_size):
for tp1, tp2 in wp[points_idx] * hop_size / fs:
    # get position on axis for a given index-pair
    coord1 = trans_figure.transform(ax1.transData.transform([tp1, 0]))
    coord2 = trans_figure.transform(ax2.transData.transform([tp2, 0]))

    # draw a line
    line = matplotlib.lines.Line2D((coord1[0], coord2[0]),
                                   (coord1[1], coord2[1]),
                                   transform=fig.transFigure,
                                   color='r')
    lines.append(line)

fig.lines = lines
plt.tight_layout()

'''
to do
1. Extract a rolling frame about the length of an average cough
2. Go frame by frame an calculate variable 'D'
3. For each frame you have a value 'D', plot ROC
4. Figure out a good threshold from ROC

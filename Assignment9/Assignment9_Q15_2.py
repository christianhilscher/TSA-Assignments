import numpy as np
import pandas as pd
from numba import njit
import matplotlib.pyplot as plt

import time
###############################################################################
@njit
def make_xmatrix(x):

    n = len(x)
    xmat = np.zeros((n,n))

    for i in range(n):

        xtmp = x[0:n-i]

        tmp = np.zeros(n)
        tmp[n-len(xtmp):] = xtmp

        xmat[i,:] = tmp
    return xmat

@njit
def make_bvector(b, x):
    n = len(x)
    bvector = np.repeat(b, n) ** np.arange(n)
    bvector = bvector * ((-1) ** np.arange(n))

    return bvector

@njit
def make_ehat_mat(b, x):
    n=len(x)

    bvec = make_bvector(b, x)
    xmat = make_xmatrix(x)

    ehat_mat = np.empty((n,n))
    for i  in np.arange(n):
        ehat_mat[:,i] = np.multiply(xmat[:,i], bvec)
    return ehat_mat

@njit
def make_ehat(ehatmat):

    n = ehatmat.shape[1]
    res = np.empty(n)
    for i in range(n):
        res[i] = np.sum(ehatmat[:,i])
    return res

@njit
def get_sigmahat(b,x):
    n = len(x)

    ehatmat = make_ehat_mat(b,x)
    eh = make_ehat(ehatmat)

    sigma = (1/(n-1)) * np.sum((eh**2))
    return sigma

@njit
def objf_MLE(b, x):
    n = len(x)
    sigmahat = get_sigmahat(b,x)
    logpart = (1/np.sqrt(2*np.pi*sigmahat))

    ehatmat = make_ehat_mat(b,x)
    eh = make_ehat(ehatmat)

    value = (n-1) * np.log(logpart) - (1/(2*sigmahat)) * np.sum((eh**2))
    return value

@njit
def calc(blist, x):
    samplesize = np.shape(x)[0]
    res_df = np.empty((samplesize, 2))

    for j in np.arange(samplesize):
        data_used = x[j,:]

        res_tmp = np.empty(len(blist))

        for i in np.arange(len(blist)):

            b_tmp = blist[i]
            res_tmp[i] = objf_MLE(b_tmp, data_used)

        b_opt = blist[np.where(np.max(res_tmp)==res_tmp)][0]
        sigma_opt = get_sigmahat(b_opt, data_used)

        res_df[j,:] = [b_opt, sigma_opt]
    return res_df

def get_biases(dataf, model):
    bias_b = (np.mean(dataf[:,0]) - b)*1000
    bias_sigma = (np.mean(dataf[:,1]) - sigma_list[model-1])*1000

    sigma_b = np.var(dataf[:,0])
    sigma_sigma = np.var(dataf[:,1])

    MSE_b = ((bias_b/1000)**2 + sigma_b) * 1000
    MSE_sigma = ((bias_sigma/1000)**2 + sigma_sigma) * 1000

    return [bias_b , bias_sigma, MSE_b, MSE_sigma]
###############################################################################

# Loading dataset - taking the one genereated in R
df_i = pd.read_csv("model_i.csv")
df_ii = pd.read_csv("model_ii.csv")
df_iii = pd.read_csv("model_iii.csv")

for df in [df_i, df_ii, df_iii]:
    # Removing index column, too lazy to do it when reading in
    df.drop('Unnamed: 0', axis=1, inplace=True)

# Converting to array
df1 = df_i.to_numpy()
df2 = df_ii.to_numpy()
df3 = df_iii.to_numpy()

b = 0.25
b_list = np.arange(-0.9,0.9, 0.05)
sigma_list = np.array([1, 5/3, 1/12]) # Providing the actual values of sigma

# Running the calculation and timing it - on my machine took roughly 1h
start_time = time.time()
# # Initialization
# c = calc(b_list, df_try)
a1 = calc(b_list, df1)
a2 = calc(b_list, df2)
a3 = calc(b_list, df3)
time.time() - start_time

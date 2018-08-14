# We import the necessary packages
from functions_global import *
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

# Cost function
costfn = rmse

# We start by loading the raw data and tidying it up.
# It contains 4 data points per vaccine configuration per time point.
dataRaw = pd.read_csv('../../data/VNA.csv')
timesData = dataRaw['days'].tolist() # List of time points
nMeasurements = 4

# We construct the arrays of data for each vaccine configuration
PBS = []  # Non-adjuvanted vaccine
MF59 = []  # Vaccine with MF59
AS03 = []  # Vaccine with AS03
Diluvac = []  #Vaccine with Diluvac
X_data = [] # List of (repeated) time points

for i in range(len(timesData)):
    for j in range(1,nMeasurements+1):
        X_data.append(timesData[i])
        PBS.append(dataRaw.T.iloc[j][i])
#    for j in range(nMeasurements+1,2*nMeasurements+1):
#        MF59.append(dataRaw.T.iloc[j][i])
#    for j in range(2*nMeasurements+1,3*nMeasurements+1):
#        AS03.append(dataRaw.T.iloc[j][i])
#    for j in range(3*nMeasurements+1,4*nMeasurements+1):
#        Diluvac.append(dataRaw.T.iloc[j][i])
PBS = np.column_stack((X_data, PBS))
#MF59 = np.column_stack((X_data, MF59))
#AS03 = np.column_stack((X_data, AS03))
#Diluvac = np.column_stack((X_data, Diluvac))

# Boundary values for the parameters to be estimated in the base case
#                gammaNA     gammaHA        mu         dmax
bounds_PBS = [(0.1, 2.5), (0.1, 7.5), (0.2, 1.0), (0.1, 0.3)]

nTotal = len(X_data)
nDays = nTotal/nMeasurements
nRuns = 2500

# indices of the datapoints separated by day
perDay = [np.arange(nMeasurements*i, nMeasurements*(i+1)) for i in range(nDays)]

def oneRun():
    # select a random sample of 1 datapoint per day
    idx_sample = [np.random.choice(x) for x in perDay]

    sample = PBS[idx_sample]
    args = (sample[:,0], sample[:,1])
    estimation = differential_evolution(costfn, bounds_PBS, args=args)
    params = estimation.x
    fun_base = estimation.fun

    return params

gammaNA_list = []
gammaHA_list = []
mu_list = []
dmax_list = []
for run in range(nRuns):
    params = oneRun()
    gammaNA, gammaHA, mu, dmax = params

    gammaNA_list.append(gammaNA)
    gammaHA_list.append(gammaHA)
    mu_list.append(mu)
    dmax_list.append(dmax)

best_fit_params = {'gammaNA': gammaNA_list, 'gammaHA': gammaHA_list,
                   'mu': mu_list, 'dmax': dmax_list}

best_fit_params = pd.DataFrame(best_fit_params)
best_fit_params.to_csv('../../params/bootstrap_base.csv', index=False)

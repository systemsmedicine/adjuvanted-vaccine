# We import the necessary packages
from functions_global import *
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

# Do we carry out the parameter fiting in log scale?
logscale = False
if logscale:
    costfn = rmse_log
    costfn_adj = rmse_adj_log
else:
    costfn = rmse
    costfn_adj = rmse_adj

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
    for j in range(nMeasurements+1,2*nMeasurements+1):
        MF59.append(dataRaw.T.iloc[j][i])
    for j in range(2*nMeasurements+1,3*nMeasurements+1):
        AS03.append(dataRaw.T.iloc[j][i])
    for j in range(3*nMeasurements+1,4*nMeasurements+1):
        Diluvac.append(dataRaw.T.iloc[j][i])
PBS = np.column_stack((X_data, PBS))
MF59 = np.column_stack((X_data, MF59))
AS03 = np.column_stack((X_data, AS03))
Diluvac = np.column_stack((X_data, Diluvac))

# Boundary values for the parameters to be estimated in the adjuvanted case
#                betaNA       betaHA       betaAb
bounds_adj = [(1.0, 5.0), (1.0, 5.0), (1.0, 60.0)]

nTotal = len(X_data)
nDays = nTotal/nMeasurements
nRuns = 2500

# indices of the datapoints separated by day
perDay = [np.arange(nMeasurements*i, nMeasurements*(i+1)) for i in range(nDays)]

# base parameters
params_base = pd.Series.from_csv('../../params/best_fit_params_base_A.csv')
gammaNA, gammaHA, mu, dmax = params_base['gammaNA'], params_base['gammaHA'], params_base['mu'], params_base['dmax']

baseQ = vQ0(np.abs(grid), dmax) + vQ0(np.abs(1 - grid), dmax)
H = Htilde*0.5*(np.sign(grid - 0.99*dmax) + np.sign(1.0 - 0.99*dmax - grid))
base_args = [gammaNA, gammaHA, mu, dmax, baseQ, H]

def oneRun():
    # select a random sample of 1 datapoint per day
    idx_sample = [np.random.choice(x) for x in perDay]
    params = np.array([])

    sample = MF59[idx_sample]
    args = (sample[:,0], sample[:,1], base_args)
    estimation = differential_evolution(costfn_adj, bounds_adj, args=args)
    params = np.append(params, estimation.x)

    sample = AS03[idx_sample]
    args = (sample[:,0], sample[:,1], base_args)
    estimation = differential_evolution(costfn_adj, bounds_adj, args=args)
    params = np.append(params, estimation.x)

    sample = Diluvac[idx_sample]
    args = (sample[:,0], sample[:,1], base_args)
    estimation = differential_evolution(costfn_adj, bounds_adj, args=args)
    params = np.append(params, estimation.x)

    return params

betaNA_M_list = []
betaHA_M_list = []
betaAb_M_list = []
betaNA_A_list = []
betaHA_A_list = []
betaAb_A_list = []
betaNA_D_list = []
betaHA_D_list = []
betaAb_D_list = []
for run in range(nRuns):
    params = oneRun()
    betaNA_M, betaHA_M, betaAb_M, betaNA_A, betaHA_A, betaAb_A, betaNA_D, betaHA_D, betaAb_D = params
    betaNA_M_list.append(betaNA_M)
    betaHA_M_list.append(betaHA_M)
    betaAb_M_list.append(betaAb_M)
    betaNA_M_list.append(betaNA_A)
    betaHA_M_list.append(betaHA_A)
    betaAb_M_list.append(betaAb_A)
    betaNA_M_list.append(betaNA_D)
    betaHA_M_list.append(betaHA_D)
    betaAb_M_list.append(betaAb_D)

best_fit_params_adj = {'betaNA_MF59': betaNA_M_list, 'betaHA_MF59': betaHA_M_list,
                       'betaAb_MF59': betaAb_M_list,
                       'betaNA_AS03': betaNA_A_list, 'betaHA_AS03': betaHA_A_list,
                       'betaAb_AS03': betaAb_A_list,
                       'betaNA_Diluvac': betaNA_D_list, 'betaHA_Diluvac': betaHA_D_list,
                       'betaAb_Diluvac': betaAb_D_list}

best_fit_params_adj = pd.DataFrame(best_fit_params_adj)
best_fit_params_adj.to_csv('../../params/bootstrap_adj.csv', index=False)

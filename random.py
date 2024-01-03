## OL-CCS
# Open-Loop Coherent Compressed Sensor
# Author: M. D. Sudeera H. Gunathilaka

##BSD 3-Clause License

#Copyright (c) 2023, NTT Research Inc., PHI labs, algorithms

#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import torch
import time
import numpy as np
import matplotlib.pyplot as plt
import math
from statistics import mean
from numpy import genfromtxt

N = 2000 # System size (N)
M = 1200 # Compressed size (N*compression ratio - This case 2000*0.6=1200)
K = 1
sparseness = 0.2 # Ratio of non-zero elemnts in the initial vector
probs = 10 # Number of problems
rmseArr = []

stt = time.time()
for pp in range(probs):
    
    torch.manual_seed(pp+123)
    AA = torch.normal(mean=0, std=1, size=(M,N), dtype=torch.double) # Observation matix
    AA = AA/math.sqrt(M)
    signal = torch.normal(mean=0, std=1, size=(N,), dtype=torch.double) # Signal (Continuous)
    support = torch.ones(N, dtype=torch.double) # Support (Binary)

    # Creating a sparse Support according to the sparseness ratio
    for i in range(N):
        if (i <= int(sparseness * N)):
            support[i] = 1.0
        else:
            support[i] = 0.0
            
    nu = 0.0
    obNoise = torch.normal(mean=0, std=0.1, size=(M,), dtype=torch.double) # Observation noise
    y = torch.matmul(AA, signal * support) + nu * obNoise # Observation signal
    A_norm = torch.zeros(N, dtype=torch.double)

    # Norm of the observation matrix
    for i in range(N):
        A_norm[i] = 0
        for j in range(M):
            A_norm[i] += AA[j][i] * AA[j][i]
        A_norm[i] = torch.sqrt(A_norm[i])

    J = torch.zeros((N,N), dtype=torch.double)
    hz = torch.zeros(N, dtype=torch.double)


    J = torch.matmul((AA/A_norm).T, AA/A_norm) # Coupling matrix
    J.fill_diagonal_(0)
    J = -J
    hz = torch.matmul((AA/A_norm).T, y) # Zeeman term

    eta_int = 0.6 # Starting L0-Regularisation parameter
    eta_end = 0.01 # Starting L0-Regularisation parameter
    REP = 50 # Number of alternating minimisation processes
    JJ = 0.25 # Feedback strength
    Tmax = 5 # Max-time for CIM
    dt = 0.1 # Time-step for CIM
    dt1 = 0.1 # Time-step for CDP
    RTmax = 10 # Max-time for CDP
    g = 1e-7 # Saturation parameter 
    Pmax = 1.5 # Max Pump rate

    sigCIM = torch.zeros(N, dtype=torch.double)
    rCDP = torch.zeros(N, dtype=torch.double)

    for ll in range(REP+1):
        eta = max(eta_int * (1 - (ll/REP)), eta_end) # L0-Regularisation parameter schedulling

        sigCIM = torch.zeros(N, dtype=torch.double)
        c = torch.zeros(N, dtype=torch.double)
        s = torch.zeros(N, dtype=torch.double)
        rr = torch.zeros(N, dtype=torch.double)

        for kk in range(round(Tmax/dt)):
            tmp = kk*dt/Tmax
            pump = tmp*tmp*Pmax # Pump rate schedulling

            # Injection field Calculation
            out = rCDP*sigCIM
            Zeeman = hz
            h_i = torch.matmul(J,out) + Zeeman

            # OL-CCS SDEs 
            w1 = torch.normal(mean=0, std=1, size=(N,), dtype=torch.double)
            w2 = torch.normal(mean=0, std=1, size=(N,), dtype=torch.double)
            rr = torch.pow(c,2) + torch.pow(s,2)
            c =  c + ((-1 + pump - rr) * c + JJ * (torch.abs(h_i) - eta)) * dt + math.sqrt(dt) * torch.sqrt(rr + 0.5) * g * w1
            s =  s + ((-1 - pump - rr) * s) * dt + math.sqrt(dt) * torch.sqrt(rr + 0.5) * g * w2
            sigCIM = ((torch.div(c, torch.abs(c))+1)/2) # Support estimation
        
        # Jacobi method for Signal optimisation
        for rr in np.arange(0.0,RTmax/dt1,1.0):

            hCDP = torch.matmul(J,(rCDP*sigCIM)) + hz
            rCDP = rCDP + (-rCDP + hCDP)*dt1

    # Evaluation by Root-Mean-Squared-Error (RMSE)
    rmse = torch.sqrt(1/N * torch.sum(torch.pow(torch.sub((sigCIM * rCDP)/A_norm, signal * support),2)))
    rmseArr.append(rmse.detach().cpu().numpy())
    print("Problem = "+ "{:.0f}".format(pp) +" | "+ " RMSE = "+ "{:.6f}".format(rmse))

print(time.time() - stt, "s")
print("==> N = "+ str(N) + ", M = "+ str(M) + " | Problems = "+ "{:.0f}".format(probs) +" | "+ " average RMSE = "+ "{:.6f}".format(sum(rmseArr)/probs))

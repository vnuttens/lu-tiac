import os
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import csv
import lmfit

# Victor Nuttens
# V02
# 24/02/2022
#--------------------------------------------------------------------------

# IN TE VULLEN

A1KL = 0.02        # AT1/A0 in left kidney
A2KL = 0.005       # AT2/A0 in left kidney
A1KR = 0.018       # AT1/A0 in right kidney
A2KR = 0.004       # AT2/A0 in right kidney

T1 = 1         # days after administration (14h d0 en 8h d1 is  0.75)
T2 = 7         # days after administration (14h d0 en 8h d7 is  6.75)
#--------------------------------------------------------------------------

a = 0
b = 8
N = 1000
t = np.linspace(a, b, N)
Teff1 = 5
Teff2 = 0.02
Teff3 = 1
Teff4 = 1

x = [0,T1,T2]
ykl = [0,A1KL,A2KL]
ykr = [0,A1KR,A2KR]

# Voor kidney, spleen, liver en speekselklier (?)
model = lmfit.models.ExpressionModel("A1*(exp(-log(2)*x/Teff1) - exp(-log(2)*x/Teff2))")
paramsKL = model.make_params(A1=A1KL, Teff1=Teff1, Teff2=Teff2)
paramsKR = model.make_params(A1=A1KR, Teff1=Teff1, Teff2=Teff2)
#params["A1"].set(vary=False)
fitkl = model.fit(ykl, paramsKL, x=x)
fitkr = model.fit(ykr, paramsKL, x=x)

Teff1kl = fitkl.params["Teff1"].value
Teff1kr = fitkr.params["Teff1"].value
Teff2kl = fitkl.params["Teff2"].value
Teff2kr = fitkr.params["Teff2"].value
A1KL = fitkl.params["A1"].value
A1KR = fitkr.params["A1"].value
print(fitkl.fit_report())
print(fitkr.fit_report())

A12kl = A1KL*(np.exp(-np.log(2)*t/Teff1kl) - np.exp(-np.log(2)*t/Teff2kl))
A12kr = A1KR*(np.exp(-np.log(2)*t/Teff1kr) - np.exp(-np.log(2)*t/Teff2kr))
TIAC_KL = np.sum(A12kl)*(b-a)/N
TIAC_KR = np.sum(A12kr)*(b-a)/N

# Voor red marrow
# To do

#--------------------------------------------------------------------------
# Making plot

fig1, ax1  = plt.subplots(1, 1, sharex=True,sharey=True)

ax1.scatter(x,ykl, color='red',label='Kidney left data')
ax1.scatter(x,ykr, color='blue',label='Kidney right data')
#ax1.plot(x,fit.init_fit, color='blue',label='Initiele curve')
ax1.plot(t,A12kl, color='red',label=f'Kidney left '
                                      f'\nTeff1 = {"{:0.2f}".format(Teff1kl)}, Teff2 = {"{:0.2f}".format(Teff2kl)}'
                                      f'\nA/A0 = {"{:0.2f}".format(TIAC_KL)}')

ax1.plot(t,A12kr, color='blue',label=f'Kidney left '
                                      f'\nTeff1 = {"{:0.2f}".format(Teff1kr)}, Teff2 = {"{:0.2f}".format(Teff2kr)}'
                                      f'\nA/A0 = {"{:0.2f}".format(TIAC_KR)}')

ax1.set_title('Time-activity curve (TAC)')
ax1.set_ylabel('A(t)/A0')
ax1.set_xlabel('Time [days]')
ax1.legend()

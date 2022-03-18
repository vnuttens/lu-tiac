import os
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import csv
import lmfit

# Victor Nuttens
# V03
# 18/03/2022
#--------------------------------------------------------------------------

# Select data file
filename = "20220309_03.csv"    # file name met data

#--------------------------------------------------------------------------

# reading csv file
rows = []
with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile,delimiter = ';')
        # extracting each data row one by one
        for row in csvreader:
                rows.append(row)


date = int(rows[1][1])
EAD = int(rows[0][1])
A0 = float(rows[8][1].replace(',', '.'))

CF = float(rows[2][1].replace(',','.'))
AT1 = int(rows[4][1].replace(',', '.'))
AT2 = int(rows[34][1].replace(',', '.'))
Vkl = float(rows[14][1].replace('Vol ','').replace(' ml','').replace(',','.'))
Vsl = float(rows[20][1].replace('Vol ','').replace(' ml','').replace(',','.'))
Vkr = float(rows[26][1].replace('Vol ','').replace(' ml','').replace(',','.'))
Vsr = float(rows[32][1].replace('Vol ','').replace(' ml','').replace(',','.'))
Csl1 = float(rows[19][1].replace('Tot ', '').replace(' CNTS','').replace(',','.'))
Csl2 = float(rows[39][1].replace('Tot ', '').replace(' CNTS','').replace(',','.'))
Csr1 = float(rows[31][1].replace('Tot ', '').replace(' CNTS','').replace(',','.'))
Csr2 = float(rows[45][1].replace('Tot ', '').replace(' CNTS','').replace(',','.'))
T1 = float(rows[3][1].replace(',','.'))
T2 = float(rows[33][1].replace(',','.'))

A1KL = Csl1/Vsl*Vkl/AT1/CF/A0
A2KL = Csl2/Vsl*Vkl/AT2/CF/A0
A1KR = Csr1/Vsr*Vkr/AT1/CF/A0
A2KR = Csr2/Vsr*Vkr/AT1/CF/A0

if not os.path.exists('Lu resultaten'):
  os.mkdir('Lu resultaten')
#--------------------------------------------------------------------------

a = 0
b = 400
N = 1000
t = np.linspace(a, b, N)
Tinit1 = 0.9
Tinit2 = 50

x = [0,T1,T2]
ykl = [0,A1KL,A2KL]
ykr = [0,A1KR,A2KR]

# Voor kidney, spleen, liver en speekselklier
model = lmfit.models.ExpressionModel("A1*(exp(-log(2)*x/Teff1) - exp(-log(2)*x/Teff2))")
paramsKL = model.make_params(A1=A1KL, Teff1=Tinit1, Teff2=Tinit2)
paramsKR = model.make_params(A1=A1KR, Teff1=Tinit1, Teff2=Tinit2)
#params["A1"].set(vary=False)
fitkl = model.fit(ykl, paramsKL, x=x)
fitkr = model.fit(ykr, paramsKR, x=x)

Teff1kl = fitkl.params["Teff1"].value
Teff1kr = fitkr.params["Teff1"].value
Teff2kl = fitkl.params["Teff2"].value
Teff2kr = fitkr.params["Teff2"].value
A1KLfit = fitkl.params["A1"].value
A1KRfit = fitkr.params["A1"].value
print(fitkl.fit_report())
print(fitkr.fit_report())

A12kl = A1KLfit*(np.exp(-np.log(2)*t/Teff1kl) - np.exp(-np.log(2)*t/Teff2kl))
A12kr = A1KRfit*(np.exp(-np.log(2)*t/Teff1kr) - np.exp(-np.log(2)*t/Teff2kr))
TIAC_KL = np.sum(A12kl)*(b-a)/N
TIAC_KR = np.sum(A12kr)*(b-a)/N

# Voor red marrow
# To do

#--------------------------------------------------------------------------
# Making plot

fig1, ax1  = plt.subplots(1, 1, sharex=True,sharey=True)

ax1.scatter(x,ykl, color='red')
ax1.scatter(x,ykr, color='blue')
#ax1.plot(x,fit.init_fit, color='blue',label='Initiele curve')
ax1.plot(t,A12kl, color='red',label=f'Kidney left '
                                      f'\nTeff1 = {"{:0.2f}".format(Teff1kl)}, Teff2 = {"{:0.2f}".format(Teff2kl)}'
                                      f'\nRes. Time = {"{:0.2f}".format(TIAC_KL)} h')

ax1.plot(t,A12kr, color='blue',label=f'Kidney left '
                                      f'\nTeff1 = {"{:0.2f}".format(Teff1kr)}, Teff2 = {"{:0.2f}".format(Teff2kr)}'
                                      f'\nRes. time = {"{:0.2f}".format(TIAC_KR)} h')

ax1.set_title(f'EAD= {EAD} \nA0 = {"{:0.1f}".format(A0)} MBq, rt = {"{:0.2f}".format(TIAC_KL+TIAC_KR)} h, mass = {"{:0.2f}".format((Vkl+Vkr)*1.04)} g')
ax1.set_xlim([-5, 200])
ax1.set_ylabel('A(t)/A0')
ax1.set_xlabel('Time [hours]')
ax1.legend()

fig1.savefig(os.path.join('Lu resultaten',f'{date}_{EAD}.png'))

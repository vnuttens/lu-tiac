import os
import matplotlib.pyplot as plt
import numpy as np
import csv
import lmfit

# Victor Nuttens
# V03
# 28/04/2022
#--------------------------------------------------------------------------
# Select data file
filename = "20220309_XX.csv"    # file name met data

#--------------------------------------------------------------------------
# Reading csv file
rows = []
with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile,delimiter = ';')
        # extracting each data row one by one
        for row in csvreader:
                rows.append(row)

EAD  = int(rows[0][1])
date2 = int(rows[1][1])
CF   = float(rows[2][1].replace(',','.'))
T3   = float(rows[3][1].replace(',','.'))
AT3  = int(rows[5][1].replace(',', '.'))
A03   = float(rows[10][1].replace(',', '.'))
Ckl3 = float(rows[17][1].replace('Tot ', '').replace(' CNTS','').replace(',','.'))
Vkl3 = float(rows[18][1].replace('Vol ','').replace(' ml','').replace(',','.'))
Ckr3 = float(rows[23][1].replace('Tot ', '').replace(' CNTS','').replace(',','.'))
Vkr3 = float(rows[24][1].replace('Vol ','').replace(' ml','').replace(',','.'))
A3KL = Ckl3/AT3/CF/A03
A3KR = Ckr3/AT3/CF/A03

#--------------------------------------------------------------------------
# Setting coefficients if this is not the first cycle

filenames1 = next(os.walk('H:\Dosimetrie Lu-177-PSMA\Lu resultaten'), (None, None, []))[2]  # [] if no file
filenames2 = next(os.walk('H:\Dosimetrie Lu-177-PSMA\Lu resultaten\in KWS'), (None, None, []))[2]  # [] if no file
filenames = np.concatenate((filenames1, filenames2))
sec_cycle = False
for file in filenames:
    if not file == 'Thumbs.db':
        if int(file[9:17]) == EAD:
            if not filename[0:8] == file[0:8]:
                if not filename[0:8] < file[0:8]:
                    filename = filename.replace(filename[0:8],file[0:8])
                    sec_cycle = True

rows = []
with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile,delimiter = ';')
        # extracting each data row one by one
        for row in csvreader:
                rows.append(row)

#--------------------------------------------------------------------------
# Setting coefficients if this is the first cycle

date = int(rows[1][1])

CF   = float(rows[2][1].replace(',','.'))
T1   = float(rows[3][1].replace(',','.'))
T2   = float(rows[4][1].replace(',','.'))
AT1  = int(rows[5][1].replace(',', '.'))
AT2  = float(rows[6][1].replace(',', '.'))
A0   = float(rows[10][1].replace(',', '.'))

Ckl1 = float(rows[17][1].replace('Tot ', '').replace(' CNTS','').replace(',','.'))
Vkl1 = float(rows[18][1].replace('Vol ','').replace(' ml','').replace(',','.'))
Ckr1 = float(rows[23][1].replace('Tot ', '').replace(' CNTS','').replace(',','.'))
Vkr1 = float(rows[24][1].replace('Vol ','').replace(' ml','').replace(',','.'))

Ckl2 = float(rows[30][1].replace('Tot ', '').replace(' CNTS','').replace(',','.'))
Vkl2 = float(rows[31][1].replace('Vol ','').replace(' ml','').replace(',','.'))
Ckr2 = float(rows[36][1].replace('Tot ', '').replace(' CNTS','').replace(',','.'))
Vkr2 = float(rows[37][1].replace('Vol ','').replace(' ml','').replace(',','.'))

A1KL = Ckl1/AT1/CF/A0
A2KL = Ckl2/AT2/CF/A0
A1KR = Ckr1/AT1/CF/A0
A2KR = Ckr2/AT1/CF/A0

if not os.path.exists('Lu resultaten'):
  os.mkdir('Lu resultaten')
#--------------------------------------------------------------------------

a = 0
b = 400
N = 1000
t = np.linspace(a, b, N)
Tinit1 = 1
Tinit2 = 1

x = [0,T1,T2]
ykl = [0,A1KL,A2KL]
ykr = [0,A1KR,A2KR]

# Voor kidney, spleen, liver en speekselklier
model = lmfit.models.ExpressionModel("A1*(exp(-log(2)*x/Teff1) - exp(-log(2)*x/Teff2))")
paramsKL = model.make_params(A1=A1KL, Teff1=Tinit1, Teff2=Tinit2)
paramsKR = model.make_params(A1=A1KR, Teff1=Tinit1, Teff2=Tinit2)

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
totalmass1 = (Vkl1+Vkr1)*1.05
totalmass2 = (Vkl2+Vkr2)*1.05
masserror = np.abs(totalmass1-totalmass2)/totalmass2
if masserror >0.2:
    totaldose = "Massa tussen tijdstip 1 en 2 verschillen meer dan 20%"
else:
    totaldose = (TIAC_KL+TIAC_KR)*A0/totalmass1*0.204*422/1000

#--------------------------------------------------------------------------
# Fitting to third point if cycle > 1

if sec_cycle:
    #errkl = A3KL/A1KL
    #errkr = A3KR/A1KR
    #x3 = [0,T3,T2]
    x3 = [0,T3]
    #ykl3 = [0,A3KL,A2KL*errkl]
    ykl3 = [0,A3KL]
    #ykr3 = [0,A3KR,A2KR*errkr]
    ykr3 = [0,A3KR]

    model = lmfit.models.ExpressionModel("A1*(exp(-log(2)*x/Teff1) - exp(-log(2)*x/Teff2))")
    paramsKL = model.make_params(A1=A1KLfit, Teff1=Teff1kl, Teff2=Teff2kl)
    paramsKR = model.make_params(A1=A1KRfit, Teff1=Teff1kr, Teff2=Teff2kr)
    #paramsKL["A1"].set(vary=False)
    paramsKL["Teff1"].set(vary=False)
    #paramsKL["Teff2"].set(vary=False)
    #paramsKR["A1"].set(vary=False)
    paramsKR["Teff1"].set(vary=False)
    #paramsKR["Teff2"].set(vary=False)
    fitkl = model.fit(ykl3, paramsKL, x=x3)
    fitkr = model.fit(ykr3, paramsKR, x=x3)

    Teff1kl3 = fitkl.params["Teff1"].value
    Teff1kr3 = fitkr.params["Teff1"].value
    Teff2kl3 = fitkl.params["Teff2"].value
    Teff2kr3 = fitkr.params["Teff2"].value
    A1KLfit3 = fitkl.params["A1"].value
    A1KRfit3 = fitkr.params["A1"].value

    A13kl = A1KLfit3 * (np.exp(-np.log(2) * t / Teff1kl3) - np.exp(-np.log(2) * t / Teff2kl3))
    A13kr = A1KRfit3 * (np.exp(-np.log(2) * t / Teff1kr3) - np.exp(-np.log(2) * t / Teff2kr3))
    TIAC_KL3 = np.sum(A13kl)*(b-a)/N
    TIAC_KR3 = np.sum(A13kr)*(b-a)/N
    totalmass3 = (Vkl3+Vkr3)*1.05
    masserror = np.abs(totalmass1-totalmass2)/totalmass2
    if masserror >0.2:
        totaldose = "Massa tussen tijdstip 1 en 3 verschillen meer dan 20%"
    else:
        totaldose = (TIAC_KL3+TIAC_KR3)*A0/totalmass3*0.204*422/1000


#--------------------------------------------------------------------------
# Making plot

fig1, ax1  = plt.subplots(1, 1, sharex=True,sharey=True,figsize=(11.69,8.27))

if sec_cycle:
    ax1.scatter(x3, ykl3, color='red')
    ax1.scatter(x3, ykr3, color='blue')
    ax1.plot(t, A13kl, color='red', label=f'Kidney left'
                                          f'\nTeff1 = {"{:0.2f}".format(Teff1kl3)}, Teff2 = {"{:0.2f}".format(Teff2kl3)}'
                                          f'\nRes. Time = {"{:0.2f}".format(TIAC_KL3)} h')
    ax1.plot(t, A13kr, color='blue', label=f'Kidney right'
                                           f'\nTeff1 = {"{:0.2f}".format(Teff1kr3)}, Teff2 = {"{:0.2f}".format(Teff2kr3)}'
                                            f'\nRes. Time = {"{:0.2f}".format(TIAC_KR3)} h')

    #ax1.scatter(x, ykl, color='red')
    #ax1.scatter(x, ykr, color='blue')
    ax1.plot(t, A12kl, color='red', linestyle='--', label=f'Kidney left, cycle 1')
    ax1.plot(t, A12kr, color='blue', linestyle='--', label=f'Kidney right, cycle 1 ')
    title = f'EAD= {EAD}, A0 = {"{:0.1f}".format(A03)} MBq \nExtrapolated kidney dose = {"{:0.2f}".format(totaldose)} Gy'
else:
    ax1.scatter(x, ykl, color='red')
    ax1.scatter(x, ykr, color='blue')

    ax1.plot(t, A12kl, color='red', label=f'Kidney left '
                                          f'\nTeff1 = {"{:0.2f}".format(Teff1kl)}, Teff2 = {"{:0.2f}".format(Teff2kl)}'
                                          f'\nRes. Time = {"{:0.2f}".format(TIAC_KL)} h')

    ax1.plot(t, A12kr, color='blue', label=f'Kidney right '
                                           f'\nTeff1 = {"{:0.2f}".format(Teff1kr)}, Teff2 = {"{:0.2f}".format(Teff2kr)}'
                                           f'\nRes. time = {"{:0.2f}".format(TIAC_KR)} h')
    title = f'EAD= {EAD}, A0 = {"{:0.1f}".format(A0)} MBq \nKidney dose = {"{:0.2f}".format(totaldose)} Gy'


ax1.set_title(title)
ax1.set_xlim([-5, 200])
ax1.set_ylabel('A(t)/A0')
ax1.set_xlabel('Time [hours]')
ax1.legend()

if sec_cycle:
    fig1.savefig(os.path.join('Lu resultaten', f'{date2}_{EAD}.pdf'))
else:
    fig1.savefig(os.path.join('Lu resultaten', f'{date}_{EAD}.pdf'))

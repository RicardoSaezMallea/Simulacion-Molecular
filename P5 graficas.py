import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import quad
import matplotlib.animation as animation

from scipy.stats import norm
from scipy.stats import linregress
from scipy.optimize import curve_fit

import matplotlib.font_manager as fm

max_it=5e6
N_ventana=int(max_it/1000)

filename = "C:/Users/ricar/Desktop/Pruebas VSC/Prubea/EnergiasVerlet_eta10_h0.01.txt"

data = np.loadtxt(filename)

T=data[:,1]
V=data[:,2]
E=data[20:,3]

#Promedios en ventanas disjuntas
T_prom=np.zeros(int(max_it/N_ventana))
V_prom=np.zeros(int(max_it/N_ventana))
E_prom=np.zeros(int(max_it/N_ventana))

for i in range(int(max_it/N_ventana)):
  T_prom[i]=np.mean(T[i*N_ventana:(i+1)*N_ventana])
  V_prom[i]=np.mean(V[i*N_ventana:(i+1)*N_ventana])
  E_prom[i]=np.mean(E[i*N_ventana:(i+1)*N_ventana])

plt.figure(figsize=(12,6))
plt.plot(T_prom,color = plt.cm.gnuplot2(0.3))
plt.plot(V_prom,color = plt.cm.gnuplot2(0.65))

plt.ylim(0.3,0.7)


plt.title("Verlet $\eta=10$ $h=0.01$",fontsize = 20)
plt.xlabel('t', fontsize = 20)
plt.ylabel('Energía', fontsize = 20)
plt.legend(("T","V"))


plt.show

V_med=np.mean(V)
T_med=np.mean(T)
#E_med=np.mean(E)

print(V_med)
print(T_med)
#print(E_med)
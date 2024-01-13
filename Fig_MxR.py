from scipy import interpolate
import math
import numpy
import pandas 

import matplotlib.pyplot as plt

from numpy import *
from pylab import *
from matplotlib import rc

# Enable usage of real TeX for labels and captions
rc('text', usetex=True)

plt.rcParams.update({
    'font.size' : 20,                   # Set font size to 11pt
    'axes.labelsize': 20,               # -> axis labels
    'legend.fontsize': 18,              # -> legends
})

def load_RNS_tab(rns_file):

    rns_table = pandas.read_csv(rns_file, sep=' ', engine='python')
  
    ratio = rns_table['ratio']
    e_15 = rns_table['e_15']
    p_15 = rns_table['p_15']
    M = rns_table['M']
    M_0 = rns_table['M_0']
    R = rns_table['r_star']
    Omega = rns_table['spin']
    Omega_k = rns_table['Omega_K']
    I = rns_table['I']
    J_M2 = rns_table['J/M^2' ]
  
    return ratio, e_15, p_15, M, M_0, R, Omega, Omega_k, I, J_M2

name_eos = ['KDE0v1'] 
name_files = ['_0', '_0.1_Omega_k', '_0.2_Omega_k', '_0.3_Omega_k', '_0.4_Omega_k', '_0.5_Omega_k', '_0.6_Omega_k', '_0.7_Omega_k', '_0.8_Omega_k', '_0.9_Omega_k', '_Max', '_2_pi_1000']

name_Omega = [r'$f = 0$', r'$f = 0.1 \: f_k$', r'$f = 0.2 \: f_K$', r'$f = 0.3 \: f_K$', r'$f = 0.4 \: f_K$', r'$f = 0.5 \: f_K$', r'$f = 0.6 \: f_K$', r'$f = 0.7 \: f_K$', r'$f = 0.8 \: f_K$', r'$f = 0.9 \: f_K$', r'$f = f_K$', r'$f = 1$ k Hz']
name_plot = [r'$KDE0v1$'] 

# M x R - EoS's - Fig 1 - a

for j in range(len(name_eos)):
        
    plt.figure(figsize=(11.0, 8.0))
    #plt.grid()
    #plt.axes().set_aspect(5.5)
    plt.xlabel(r'$R_e$ (km)') 
    plt.ylabel(r'$M (M_\odot)$')
    
    M_points = {}
    R_points = {}
    
    for i in range(len(name_files)):
        
        eos_name = "C:/Users/carol/Analise_de_dados_py/Output_Fig_MxR/" + name_eos[j] + name_files[i] + ".txt"
        
        ratio, e_15, p_15, M, M_0, R, Omega, Omega_k, I, J_M2 = load_RNS_tab(eos_name)
        
        y = p_15/e_15
        
        R = numpy.array(R)
        M = numpy.array(M)
        y = numpy.array(y)
        
        m = len(M[M <= 1])
        
        M = M[m:]
        R = R[m:]
        y = y[m:]
        
        for n in range(1, len(M)):
            
            dif = M[n] - M[n-1] 
            
            if dif < 0:
                
                M = M[:n+1]
                R = R[:n+1]
                y = y[:n+1]
                
                break
        
        R = numpy.flip(R)
        M = numpy.flip(M)
        y = numpy.flip(y)
        
        M_x_R = interpolate.interp1d(R, M)
        R_values = numpy.arange(R[0], R[-1], 0.01)
        
        M_values = [M_x_R(R_values[n]) for n in range(0, len(R_values))]
        
        # Curva f = 1000 Hz #
        
        if name_files[i] == '_2_pi_1000':
            
            plt.plot(R_values, M_values, ':', color = 'grey', label = name_Omega[i])
            
        # Demais curvas #
            
        else:
            
            # P_c/rho_c = cte curves #
            
            R_x_y = interpolate.interp1d(y, R)
            M_x_y = interpolate.interp1d(y, M)
            
            #R_critical = R_x_y(0.3)
            #M_critical = M_x_y(0.3)
            
            #C_critical = M_critical/R_critical # Keep the dimension for this plot
            
            y_values = numpy.arange(0.1, 0.40, 0.03)
            y_values = numpy.flip(y_values)
            
            R_values_2 = [R_x_y(y_values[u]) for u in range(0,len(y_values))]
            M_values_2 = [M_x_y(y_values[u]) for u in range(0,len(y_values))]
            
            R_points["R" + str(name_eos[j]) + str(i)] = R_values_2
            M_points["M" + str(name_eos[j]) + str(i)] = M_values_2
            
            # Curva para Omega = 0 #
            
            if name_files[i] == '_0':
            
                plt.plot(R_values, M_values, color = 'black', label = name_Omega[i]) 
            
            # Curva para Omega = 0.9 Omega_K #
            
            if name_files[i] == '_0.9_Omega_k':
                
                plt.plot(R_values, M_values, color = 'gray', label = name_Omega[i])
            
            # Curva para Omega = Omega_K #
                
            if name_files[i] == '_Max':
                
                #C_Max_critical = C_critical
                #R_Max_critical = R_critical
                
                M_Max = M[0] + 0.05 # Define o limite superior do eixo y no plot
                    
                plt.plot(R_values, M_values, '--', color = 'black', label = r'$f = f_k$')
    
    for g in range(len(name_Omega)-1):
    
        for u in range(len(M_points["M" + str(name_eos[j]) + str(g)])):
            
            M_plot_points = numpy.zeros(len(name_Omega)-1)
            R_plot_points = numpy.zeros(len(name_Omega)-1)
     
            for v in range(len(name_Omega)-1):
                
                M_plot_points[v] = M_points["M" + str(name_eos[j]) + str(v)][u]
                R_plot_points[v] = R_points["R" + str(name_eos[j]) + str(v)][u]
            
            # Curvas p_c/rho_c = cte #
            
            New_M_x_R = interpolate.interp1d(R_plot_points, M_plot_points)
            
            Final_R = numpy.arange(R_plot_points[0], R_plot_points[-1], 0.1)
            Final_M = [New_M_x_R(Final_R[n]) for n in range(0,len(Final_R))]

            plt.plot(R_plot_points, M_plot_points, '--', color = 'darkgray')
            
            # Retas para Omega = 0 #
            
            x = [0, R_plot_points[0]]
            y = [0, M_plot_points[0]]
            
            x_0 = R_points["R" + str(name_eos[j]) + str(0)][0] - 0.8
            x_f = R_points["R" + str(name_eos[j]) + str(len(name_Omega)-2)][-1] + 0.8
            
            yxx = interpolate.interp1d(x, y, fill_value = 'extrapolate')
            
            x_values = numpy.arange(0, R_plot_points[-1] + 5, 0.1)
            y_values = [yxx(x_values[n]) for n in range(0,len(x_values))]
            
            plt.plot(x_values, y_values, linewidth = 0.1, color = 'lightgray')
    
    plt.xlim(x_0 - 1, x_f + 1) # Define o limite para o eixo x 
    plt.ylim(1.4, M_Max + 0.15) # Define o limite para o eixo y

    plt.legend(loc = 2)
    plt.text((x_0 + x_f)/2, M_Max + 0.07, name_plot[j], ha='center', va='center', fontsize=20)
    plt.show()
    
    
    
    
    
    
    
#*- coding: utf-8 -*-
# .. V1.0 December 11, 2025 - Massimo Donelli - Routine for read the Nano VNA ..
# .. V1.0 July 21, 2026 - Massimo Donelli - Filter added to remove reading noise ..
# .. Project for Dielectric Characteriaztion with an open end coaxial probe ..
# -- Reference Paper:
# "Dielectric Permettivity Measurement Using Open-Ended Coaxial Probe - Modeling and Simulation Based on the SImple Capacitive Modelproperties of minerals at microwave heating frequencies using an open-ended 
# coaxial line," by Antonio Saloric, and Andela Matkovic, Sensor DOI 10.3390/s22166024
# .. The fabricated probe consists of a alluminium external pipe dim 12.4mm and an internal brass pipe of dim. 4mm the 
# .. theoretical values of Z0 = 67,88 Ohms, the measured value with VNA results 60,11 Ohms
# .. Copyright@2025
# .. University of Trento, 38100 Trento, Italy
import pynanovna
import matplotlib.pyplot as plt
import math
from scipy.signal import savgol_filter
import numpy as np
import os
#
def INITIALIZE(BEGIN_FREQ,END_FREQ,POINTS_NUMBER):
    VNA = pynanovna.VNA()
    VNA.set_sweep(BEGIN_FREQ, END_FREQ, POINTS_NUMBER)
    return VNA
# .. Read the S parameters from NanoVNA ..
def READ_S_PARAMETERS(DEVICE, FLAG_PLOT):
    S11, S21, FREQ = DEVICE.sweep()
    if(FLAG_PLOT == 1):
        # .. Create a figure panel ..
        fig, axs = plt.subplots(2)
        # .. Set the figure axis ..
        MIN_Y = min(min(S11.real),min(S11.imag))
        MAX_Y = max(max(S11.real),max(S11.imag))
        axs[0].set_xlim(min(FREQ),max(FREQ))
        axs[0].set_ylim(MIN_Y,MAX_Y+MAX_Y*0.1)
        axs[0].plot(FREQ, S11.real,'-',linewidth='1.5',color='blue')
        axs[0].plot(FREQ, S11.imag,'-',linewidth='1.5',color='red')
        axs[0].set_title('$S_{11}$')
        axs[0].legend(['$Re[S_{11}]$', '$Im[S_{11}]$'])
        axs[0].set_xlabel("Frequency")
        axs[0].grid()
        # .. Now plot the |S11| .. #
        MIN_Y = min(abs(S11))
        MAX_Y = max(abs(S11))
        axs[1].set_xlim(min(FREQ),max(FREQ))
        axs[1].set_ylim(MIN_Y,MAX_Y+MAX_Y*0.1)
        axs[1].plot(FREQ, abs(S11),'-',linewidth='1.5',color='red')
        axs[1].set_title('$|S_{11}|$')
        axs[1].legend(['$|S_{11}|$'])
        axs[1].set_xlabel("Frequency")
        axs[1].set_ylabel('$|S_{11}|$ [dB]')
        axs[1].grid()
        plt.show()
    return FREQ, S11, S21
# .. This function visualize the interpolated data ..
def VisualizeResults(X,Y,Z):
    # .. Plot Synthetic Data ..    
        plt.plot(X, Y,'-',linewidth='1.5',color='red')
        plt.plot(X, Z,'-',linewidth='1.5',color='blue')
        plt.grid()  
        plt.title("Predicted Data")
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Dielectric permittivity')
        plt.legend(['Re(EPS)','Im(EPS)'])
        plt.show()
# .. This function generate the data to be interpolated by the SVM ..
def ReadData(FILE_NAME_RE_EPS,FILE_NAME_IM_EPS,FLAG_VISUALIZATION):
    DUMMY_FREQ_RE = []
    DUMMY_FREQ_IM = []
    DUMMY_RE_EPS = []
    DUMMY_IM_EPS = []
    # .. Open the files ..
    FILE_RE_EPS = open(FILE_NAME_RE_EPS,'r')
    FILE_IM_EPS = open(FILE_NAME_IM_EPS,'r')
    # .. Read the real part of EPS Water ..
    for LINEE in FILE_RE_EPS: 
        DUMMY_FREQ_RE.append(float(LINEE.split(',')[0]))
        DUMMY_RE_EPS.append(float(LINEE.split(',')[1]))
    # .. Read the real part of EPS Water ..
    for LINEE in FILE_IM_EPS:
        DUMMY_FREQ_IM.append(float(LINEE.split(',')[0]))
        DUMMY_IM_EPS.append(float(LINEE.split(',')[1]))    
    # .. Close Files
    FILE_RE_EPS.close()
    FILE_IM_EPS.close()
    # .. Convert into numpy arry ..
    FREQ_RE = np.zeros((len(DUMMY_FREQ_RE),1))
    RE_EPS = np.zeros((len(DUMMY_FREQ_RE),1))
    for I_X_FOR in range(len(DUMMY_FREQ_RE)):
        FREQ_RE[I_X_FOR] = DUMMY_FREQ_RE[I_X_FOR]
        RE_EPS[I_X_FOR] = DUMMY_RE_EPS[I_X_FOR]
    #
    FREQ_IM = np.zeros((len(DUMMY_FREQ_IM),1))
    IM_EPS = np.zeros((len(DUMMY_FREQ_IM),1))
    for I_X_FOR in range(len(DUMMY_FREQ_IM)):
        FREQ_IM[I_X_FOR] = DUMMY_FREQ_IM[I_X_FOR]
        IM_EPS[I_X_FOR] = DUMMY_IM_EPS[I_X_FOR]
    #  
    if(FLAG_VISUALIZATION==1):
        # .. Plot Synthetic Data ..    
        plt.plot(FREQ_RE, RE_EPS,'o-',linewidth='1.5',color='red')
        plt.plot(FREQ_IM, IM_EPS,'o-',linewidth='1.5',color='blue')
        #  
        plt.grid()  
        plt.title('Water Dielectric permettivity 25 °C')
        plt.xlabel('Frequency [GHz]')
        plt.legend(['Re(EPS)','Im(EPS)'])
        plt.show()     
    return FREQ_RE,RE_EPS,FREQ_IM,IM_EPS
# .. This function is mandatory to guarantee a good SVM regression ..
def ConvertFreqInt(FREQ_INT):
    LENGTH = len(FREQ_INT)
    FREQ = np.zeros((LENGTH,1))
    for I_X_FOR in range(LENGTH):
        FREQ[I_X_FOR] = FREQ_INT[I_X_FOR]*10^9
    return FREQ
# .. Retrieve the dielectric permettivity of the water at a given temperature with SVM ..
def GenerateWaterEpsValues1(FREQ,FLAG_VISUALIZATION):
    EPS_INF = 5.2
    EPS_W = 78.63
    TAU_W = 8.27*10e-12
    WATER_EPS = []
    RE_EPS = []
    IM_EPS = []
    for J_X_FOR in range(len(FREQ)):
        OMEGA = 2*math.pi*FREQ[J_X_FOR]  
        DUMMY_EPS = EPS_INF + ((EPS_W-EPS_INF)/complex(1,OMEGA*TAU_W))
        RE_EPS.append(DUMMY_EPS.real)
        IM_EPS.append(DUMMY_EPS.imag)
        WATER_EPS.append(DUMMY_EPS)   
    if(FLAG_VISUALIZATION==1):
        VisualizeResults(FREQ,RE_EPS,IM_EPS)
    return WATER_EPS   
# .. This Subroutine calibrate the fringing capacities of the coaxial probe ..
def PROBE_CALIBRATION(DEVICE,EPS_WATER,FLAG_PLOT):
    input("Connect the probe to VNA connect the open end termination, then push a button...")
    print("Read the S11 of the Open End Termination")
    FREQ, S11_OPEN, S21 = READ_S_PARAMETERS(DEVICE, FLAG_PLOT)
    input("Connect the short to the probe, then and push a bottom")
    print("Read the S11 of the Short End Termination...")
    FREQ, S11_SHORT, S21 = READ_S_PARAMETERS(DEVICE, FLAG_PLOT)
    input("Connect a water sample at known temperture to the probe, then and push a bottom")
    print("Read the S11 of the water sample...")
    # .. Initialize the NanoVNA ..
    FREQ, S11_WATER, S21 = READ_S_PARAMETERS(DEVICE, FLAG_PLOT)
    LENGTH = len(FREQ)
    A1 = []
    A2 = []
    for J_X_FOR in range(LENGTH):
        DENOM = S11_WATER[J_X_FOR] - S11_OPEN[J_X_FOR]
        A1_NUM = (EPS_WATER[J_X_FOR]*(S11_SHORT[J_X_FOR]-S11_WATER[J_X_FOR])-(S11_SHORT[J_X_FOR]-S11_OPEN[J_X_FOR]))
        A2_NUM = (EPS_WATER[J_X_FOR]*(S11_SHORT[J_X_FOR]-S11_WATER[J_X_FOR])*S11_OPEN[J_X_FOR])-S11_WATER[J_X_FOR]*(S11_SHORT[J_X_FOR]-S11_OPEN[J_X_FOR])
        A1.append(A1_NUM/DENOM)
        A2.append(A2_NUM/DENOM)
    return A1,A2,S11_SHORT # .. A3 correspond to short circuit reflection coefficient ..
# .. Save the calibration data of the probe ..
def SAVE_CALIBRATION(FREQ, A1, A2, A3, FILE_NAME_CALIBRATION):
     # .. Open File
    file = open(FILE_NAME_CALIBRATION,'w')
    # .. Write Header file
    file.write("# Freq. A1  A2  A3\n")
    # .. The abs() function has been inserted to remove the negative values of capacitance ..
    for I_X_FOR in range(len(FREQ)):
        stringa=str(FREQ[I_X_FOR])+","+str(A1[I_X_FOR])+","+str(A2[I_X_FOR])+","+str(A3[I_X_FOR])+"\n"
        file.write(stringa)
    # .. Close File
    file.close()
def READ_CALIBRATION(FILE_NAME_CALIBRATION):
    file = open(FILE_NAME_CALIBRATION,'r')
    # .. Read the real part of EPS Water ..
    A1 = []
    A2 = []
    A3 = []
    FREQ = []
    for LINE in file: 
        if(LINE.startswith("#")==0):
            FREQ.append(float(LINE.split(',')[0]))
            A1.append(complex(LINE.split(',')[1]))
            A2.append(complex(LINE.split(',')[2]))
            A3.append(complex(LINE.split(',')[3]))
    # .. Close File
    file.close()  
    return FREQ,A1,A2,A3 
def LIMIT_FREQ(FREQ,RE_EPS,IM_EPS):
    MIN_FREQ = 400e6
    DUMMY_FREQ = []
    DUMMY_RE_EPS = []
    DUMMY_IM_EPS = []
    for I_X_FOR in range(len(FREQ)):
        if(FREQ[I_X_FOR]>MIN_FREQ):
         DUMMY_FREQ.append(FREQ[I_X_FOR])
         DUMMY_RE_EPS.append(RE_EPS[I_X_FOR])
         DUMMY_IM_EPS.append(IM_EPS[I_X_FOR])
    return DUMMY_FREQ,DUMMY_RE_EPS,DUMMY_IM_EPS
#
def MEASURE_EPS(DEVICE, FILE_NAME_CALIBRATION, FLAG_PLOT):
    FREQ, A1, A2, A3 = READ_CALIBRATION(FILE_NAME_CALIBRATION)
    FREQ, S11_MUT, S21_MUT = READ_S_PARAMETERS(DEVICE, 0)
    EPS_MUT = []
    for I_X_FOR in range(len(FREQ)):
        EPS_MUT.append((A1[I_X_FOR]*S11_MUT[I_X_FOR] - A2[I_X_FOR])/(A3[I_X_FOR]-S11_MUT[I_X_FOR]))
    # .. Convert into real and imaginary EPS ..
    MUT_RE_EPS = []
    MUT_IM_EPS = []
    for I_X_FOR in range(len(FREQ)):
        MUT_RE_EPS.append(EPS_MUT[I_X_FOR].real)
        MUT_IM_EPS.append(EPS_MUT[I_X_FOR].imag)
    FREQ,MUT_RE_EPS,MUT_IM_EPS = LIMIT_FREQ(FREQ,MUT_RE_EPS,MUT_IM_EPS)
    if(FLAG_PLOT == 1):
        # .. Create a figure panel ..
        fig, axs = plt.subplots(2,2)
        # .. Now plot the EPS1 .. #
        # .. Set the figure axis ..
        MIN_Y = min(MUT_RE_EPS)
        MAX_Y = max(MUT_RE_EPS)
        axs[0][0].set_xlim(min(FREQ),max(FREQ))
        axs[0][0].set_ylim(MIN_Y,MAX_Y)
        axs[0][0].plot(FREQ, MUT_RE_EPS,'-',linewidth='1.5',color='blue')
        axs[0][0].plot(FREQ,savgol_filter(MUT_RE_EPS,20,1),'--',color='cyan')
        axs[0][0].set_ylabel("ε'")
        axs[0][0].set_xlabel("Frequency [Hz]")
        axs[0][0].legend(["Measured ε'","Measured and filtered ε'"])
        axs[0][0].grid()
        # .. Now plot the EPS2 .. #
        MIN_Y = min(MUT_IM_EPS)
        MAX_Y = max(MUT_IM_EPS)
        axs[1][0].set_xlim(min(FREQ),max(FREQ))
        axs[1][0].set_ylim(MIN_Y,MAX_Y+MAX_Y*0.1)
        axs[1][0].plot(FREQ, MUT_IM_EPS,'-',linewidth='1.5',color='red')
        axs[1][0].plot(FREQ,savgol_filter(MUT_IM_EPS,20,1),'--',color='green')
        axs[1][0].set_xlabel("Frequency [Hz]")
        axs[1][0].set_ylabel("ε\"")
        axs[1][0].legend(["Measured ε\"","Measured and filtered ε\""])
        axs[1][0].grid()
        #
        MIN_Y = min(MUT_RE_EPS)
        MAX_Y = max(MUT_RE_EPS)
        axs[0][1].set_xlim(min(FREQ),max(FREQ))
        axs[0][1].set_ylim(MIN_Y,MAX_Y)
        axs[0][1].plot(FREQ,savgol_filter(MUT_RE_EPS,20,1),'-',linewidth='1.5',color='blue')
        axs[0][1].set_ylabel("ε'")
        axs[0][1].set_xlabel("Frequency [Hz]")
        axs[0][1].legend(["Measured ε'","Measured and filtered ε'"])
        axs[0][1].grid()
        #
        MIN_Y = min(MUT_IM_EPS)
        MAX_Y = max(MUT_IM_EPS)
        axs[1][1].set_xlim(min(FREQ),max(FREQ))
        axs[1][1].set_ylim(MIN_Y,MAX_Y)
        axs[1][1].plot(FREQ,savgol_filter(MUT_IM_EPS,20,1),'-',color='red')
        axs[1][1].set_ylabel("ε'")
        axs[1][1].set_xlabel("Frequency [Hz]")
        axs[1][1].legend(["Measured ε\""])
        axs[1][1].grid()
        plt.show()
    return FREQ,EPS_MUT,savgol_filter(MUT_RE_EPS,20,1),savgol_filter(MUT_IM_EPS,20,1)
#
def SAVE_MEASURE_EPS(FILE_NAME, FREQ, EPS, F_RE_EPS, F_IM_EPS):
    file = open(FILE_NAME,'w')
    for I_X_FOR in range(len(FREQ)):
        stringa=str(FREQ[I_X_FOR])+"\t"+str(EPS[I_X_FOR].real)+"\t"+str(EPS[I_X_FOR].imag)+"\t"+str(F_RE_EPS[I_X_FOR])+"\t"+str(F_IM_EPS[I_X_FOR])+"\n"
        file.write(stringa)
    # .. Close File
    file.close()  
    # .. Read the real part of EPS Water ..
#
if __name__=="__main__":
    BEGIN_FREQ = 0.50e9
    END_FREQ = 1.4e9
    POINTS_NUMBER = 101
    FLAG_PLOT = 1
    FILE_NAME_CALIBRATION = "Calibration_M.txt"
    Z0 = 50.0 # .. Our probe Z0 measured with VNA ..
    # .. Initialize the NanoVNA ..
    DEVICE = INITIALIZE(BEGIN_FREQ,END_FREQ,POINTS_NUMBER)
    os.system('clear all')
    #
    print("_______________________________________________________________________")
    print("Do you want to calibrate C0 and Cf or to measure Eps material?")
    Flag = input("Insert C to calibrate or M to measure:")
    if((Flag == 'C') or (Flag == 'c')):
      print("Probe Calibration")
      FREQ, DUMMY_S11, DUMMY_S21 = READ_S_PARAMETERS(DEVICE, 0)
      EPS_WATER = GenerateWaterEpsValues1(FREQ,FLAG_PLOT)
      A1, A2, A3 =PROBE_CALIBRATION(DEVICE,EPS_WATER,FLAG_PLOT)
      SAVE_CALIBRATION(FREQ, A1, A2, A3, FILE_NAME_CALIBRATION)
      print("Calibration Data Saved ..")
      print("End ...")
      print("__________________________________")
    elif((Flag == 'M') or (Flag == 'm')):
        print("Eps Material Estimation")
        print("Material EPS Measurement")
        FREQ, EPS_MUT, F_RE_EPS, F_IM_EPS= MEASURE_EPS(DEVICE, FILE_NAME_CALIBRATION, FLAG_PLOT)
        Flag = input("Do you want to save the measures (Y/N):")
        if((Flag == 'Y') or (Flag == 'y')):
            FILE_NAME_MEASURE = input("Insert file name of the measure:")
            SAVE_MEASURE_EPS(FILE_NAME_MEASURE, FREQ, EPS_MUT, F_RE_EPS, F_IM_EPS)
        print("End ...")
        print("__________________________________")
    else:
        print("You inserted the Wrong Value")
        exit(0)
        print("End ...")
        print("__________________________________")
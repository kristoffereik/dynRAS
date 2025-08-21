import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ChemODE_BIO import chemODE_BIO #simulation case 1 modify uncomment function for the right hand side of for the ODE of the biofilter class to run it
from ChemODE_BIO import chemODE_BIO_HCO3 #simulation case 1 HCO3 dosing only uncomment the function for the right hand side of for the ODE of the biofilter class to run it
from ChemODE_BIO import chemODE_BIO_NaOH #simulation case 1 NaOH dosing only uncomment the function for the right hand side of for the ODE of the biofilter class to run it
from ChemODE_BIO import chemODE_BIO_alk_control
from ChemODE_BIO import chemODE_BIO_pH_control #set a pH threshold and dose NaOH or HCO3 to maintain the pH at the setpoint uncomment the function for the right hand side of for the ODE of the biofilter class to run it
from ChemODE_Fish import chemODE_FISH
from scipy.integrate import solve_ivp
from params import params
from Fish_growth import Weight
from Biomass_function import Biomass
from ChemODE_DGS import chemODE_DGS

class Sump:
    def __init__(self, CO2aq, HCO3, CO32, H, OH,NH4, NH3, NO2,V,exchange_rate):
        self.CO2aq = CO2aq
        self.HCO3 = HCO3
        self.CO32 = CO32
        self.H = H
        self.OH = OH
        self.NH4 = NH4
        self.NH3 = NH3
        self.V = V
        self.NO2 = NO2
        self.exchange_rate = exchange_rate


class Dosing_pump:
    def __init__(self, OH, dosing_compartment,V):
        self.OH = OH
        self.dosing_compartment = dosing_compartment
        self.V=V

class Dosing_pump_HCO3:
    def __init__(self, HCO3, dosing_compartment,V):
        self.HCO3 = HCO3
        self.dosing_compartment = dosing_compartment
        self.V=V


class fish_tank:
    Fish_Biomass_Max = 50000

    def __init__(self, CO2aq, HCO3, CO32, H, OH, NH4, NH3,NO2, Tminus, V,
                 Fish_weight, Fish_number, Fish_Biomass):
        self.CO2aq = CO2aq
        self.HCO3 = HCO3
        self.CO32 = CO32
        self.H = H
        self.OH = OH
        self.NH4 = NH4
        self.NH3 = NH3
        self.NO2 = NO2
        self.Tminus = Tminus
        self.V = V
        self.Fish_weight = Fish_weight
        self.Fish_number = Fish_number
        self.Fish_Biomass = Fish_Biomass
        self.CO2_fish=[]
        self.CO2_fish_t=[]

    def RHS(self, params, t):
        dY1 = chemODE_FISH(self, params,t)
        dY2 = Weight(self, params, t)
        dY = np.append(dY1, dY2)
        return dY



class biofilter:
    def __init__(self,CO2aq, HCO3, CO32, H, OH, NH4, NH3, NO2, Tminus,dosing_OH,dosing_HCO3, V, AOB,
                 NOB,B):
        self.CO2aq = CO2aq
        self.HCO3 = HCO3
        self.CO32 = CO32
        self.H = H
        self.OH =OH
        self.NH4 = NH4
        self.NH3 = NH3
        self.NO2 = NO2
        self.Tminus = Tminus
        self.dosing_OH = dosing_OH
        self.dosing_HCO3 = dosing_HCO3
        self.V = V
        self.AOB = AOB
        self.NOB = NOB
        self.B=B
        self.dosing_amount_OH=[] # store the dosing amount of OH
        self.dosing_amount_HCO3=[] # store the dosing amount of HCO3
        self.dosing_time=[]
    def RHS(self, params, t):
        dY = chemODE_BIO(self, params,t) #uncomment to run and comment the three other lines (Jafai et al. 2024 simulation settings)
        #dY = chemODE_BIO_HCO3(self, params, t) #uncomment to run and comment the three other lines (simulation case 1 -HCO3 dosing only)
        #dY = chemODE_BIO_NaOH(self, params, t) #uncomment to run and comment the three other lines (simulation case 2 -NaOH dosing only)
        #dY= chemODE_BIO_alk_control(self,params, t) #uncomment to run and comment the three other (simulation case 3 -NaOH & HCO3 dosing based on CO2)
        #dY=chemODE_BIO_pH_control(self,params, t) #uncomment to run and comment the three other lines set a threshold for the pH and control it using eiter OH or HCO3
        return dY


class degasser:
    def __init__(self, CO2aq, HCO3, CO32, H, OH, NH4, NH3, NO2, Tminus, V,Tplus):
        self.CO2aq = CO2aq
        self.HCO3 = HCO3
        self.CO32 = CO32
        self.H = H
        self.OH = OH
        self.NH4 = NH4
        self.NH3 = NH3
        self.Tminus = Tminus
        self.V = V
        self.NO2 = NO2
        self.Tplus=Tplus

    def RHS(self, params, t):
        dY1 = chemODE_DGS(self, params,t)
        return dY1

def chem(t, S, params, FishTank, B1, DGS):
    # Assigning values from the array S to the  FishTank
    FishTank.CO2aq = S[0]
    FishTank.HCO3 = S[1]
    FishTank.CO32 = S[2]
    FishTank.H = S[3]
    FishTank.OH = S[4]
    FishTank.NH4 = S[5]
    FishTank.NH3 = S[6]
    FishTank.NO2 = S[7]
    FishTank.Fish_weight = S[8]
    FishTank.Fish_Biomass = S[9]

    # Assigning values to B1 (biofilter compartment)
    B1.CO2aq = S[10]
    B1.HCO3 = S[11]
    B1.CO32 = S[12]
    B1.H = S[13]
    B1.OH = S[14]
    B1.NH4 = S[15]
    B1.NH3 = S[16]
    B1.NO2 = S[17]
    B1.AOB = S[18]
    B1.NOB = S[19]

    # Assigning values to DGS (degasser compartment)
    DGS.CO2aq = S[20]
    DGS.HCO3 = S[21]
    DGS.CO32 = S[22]
    DGS.H = S[23]
    DGS.OH = S[24]
    DGS.NH4 = S[25]
    DGS.NH3 = S[26]
    DGS.NO2 = S[27]

    # Compute RHS for each component
    dY1 = FishTank.RHS(params, t)
    dY2 = B1.RHS(params, t)
    dY3 = DGS.RHS(params, t)
    dY4 = np.append(dY1, dY2)
    dY = np.append(dY4, dY3)

    return dY
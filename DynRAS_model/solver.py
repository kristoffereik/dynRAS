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
# from params import params
from Fish_growth import Weight
from Biomass_function import Biomass
from ChemODE_DGS import chemODE_DGS

from class_definition import Sump, Dosing_pump, Dosing_pump_HCO3, fish_tank, biofilter, degasser, chem

class Solver:

    def __init__(self, params):
        self.params = params

        # when choosing the initial conditions, make sure that the initial conditions are consistent with the chemical equilibrium e.g. you cannot have a higher NH3 concentration than NH4 at pH 7
        # When running scenario 1,2,3 change the initial conditions of HCO3 to start with the targeted alkalinity to obtain the same results as in the paper
        # initial condition Fish tank
        self.CO2aq_FT = 10/44.1 #10 mg/L / 44.1 g/mol 10/44.1 mmol/L
        self.HCO3_FT =200/61.01 #(mg/l)==> mmol.l-1
        #self.HCO3_FT =70/61.01
        self.CO32_FT =10**-6*10**3
        self.H_FT = 10**-7.6*10**3
        self.OH_FT = 10**-6.4*10**3
        self.NH4_FT =0
        self.NH3_FT =0
        self.NO2_FT =0

        self.Fishweight_FT = self.params.w0

        self.Fish_number_FT = 525

        self.Fish_Biomass_FT = self.Fish_number_FT * self.Fishweight_FT

        # initial condition B1


        self.CO2aq_B1 =10/44.1
        self.HCO3_B1 =200/61.01
        #self.HCO3_B1 =70/61.01
        self.CO32_B1 =10**-6*10**3
        self.H_B1 =10**-7.6*10**3
        self.OH_B1 = 10**-6.5*10**3
        self.NH4_B1 =0
        self.NH3_B1 =0
        self.NO2_B1 =0
        self.AOB_B1 = 960000
        self.NOB_B1 = 480000
        self.B=0.3*80000*10**3*0.05
        # initial condition DGS
        self.CO2aq_DGS =10/44.1
        self.HCO3_DGS =200/61.01
        #self.HCO3_DGS =70/61.01
        self.CO32_DGS=10**-6*10**3
        self.H_DGS =10**-7.6*10**3
        self.OH_DGS =10**-6.4*10**3
        self.NH4_DGS =0
        self.NH3_DGS =0
        self.NO2_DGS =0
        self.NO3_DGS =40/62


        # initial values for the dosing pump
        #pump_OH
        self.OH_dosing=2000 #mmol/l (concentration of the dosing NaOH solution)
        #pump_HCO3
        self.HCO3_dosing=2000 #mmol/l (concentration of the dosing NaHCO3 solution)
        #value for the sump
        self.CO2aq_Sump =0
        self.HCO3_Sump =70/61.01
        self.CO32_Sump=10**-6*10**3
        self.H_Sump =10**-7.5*10**3
        self.OH_Sump =10**-6.5*10**3
        self.NH4_Sump =0
        self.NH3_Sump =0
        self.NO2_Sump =0
        self.V_Sump = 50

        self.FishTank = fish_tank(self.CO2aq_FT, self.HCO3_FT, self.CO32_FT, self.H_FT, self.OH_FT,self.NH4_FT, self.NH3_FT,
                             self.NO2_FT, [], 1000, self.Fishweight_FT, self.Fish_number_FT, self.Fish_Biomass_FT)
        self.B1 = biofilter( self.CO2aq_B1, self.HCO3_B1, self.CO32_B1, self.H_B1, self.OH_B1, self.NH4_B1, self.NH3_B1,
                     self.NO2_B1, self.FishTank,[],[], 800, self.AOB_B1, self.NOB_B1, self.B)
        self.DGS = degasser( self.CO2aq_DGS, self.HCO3_DGS, self.CO32_DGS, self.H_DGS, self.OH_DGS, self.NH4_DGS, self.NH3_DGS,
                       self.NO2_DGS, self.B1, 700,[])

        self.exchange_rate = (self.B1.V+self.DGS.V+self.FishTank.V)*0.25/(60*60*24) # 25% of the total volume is exchanged every day this can be change to any % of the total volume

        self.Sump_1= Sump( self.CO2aq_Sump, self.HCO3_Sump, self.CO32_Sump, self.H_Sump, self.OH_Sump, self.NH4_Sump, self.NH3_Sump, self.NO2_Sump,self.V_Sump, self.exchange_rate)
        self.DGS = degasser( self.CO2aq_DGS, self.HCO3_DGS, self.CO32_DGS, self.H_DGS, self.OH_DGS, self.NH4_DGS, self.NH3_DGS, self.NO2_DGS, self.B1, 700,self.Sump_1)
        self.dosing_pump_OH=Dosing_pump(self.OH_dosing,self.B1,1000)
        self.dosing_pump_HCO3=Dosing_pump_HCO3(self.HCO3_dosing,self.B1,1000)
        self.B1.dosing_OH=self.dosing_pump_OH
        self.B1.dosing_HCO3=self.dosing_pump_HCO3
        self.FishTank.Tminus = self.DGS

        self.S0 = [ self.CO2aq_FT, self.HCO3_FT, self.CO32_FT, self.H_FT, self.OH_FT, self.NH4_FT, self.NH3_FT,self.NO2_FT, self.Fishweight_FT,self.Fish_Biomass_FT, self.CO2aq_B1, self.HCO3_B1, self.CO32_B1, self.H_B1, self.OH_B1, self.NH4_B1, self.NH3_B1, self.NO2_B1, self.AOB_B1, self.NOB_B1, self.CO2aq_DGS, self.HCO3_DGS, self.CO32_DGS, self.H_DGS, self.OH_DGS,
             self.NH4_DGS, self.NH3_DGS,self.NO2_DGS]

    def solve(self):
        for cycle in range(0,10): #  solve the system for 10 cycles of 14 days, after each cycle the fish biomass is reduced back to 50kg/m3, the simulation is in seconds

            tspan = [14 * 24 * 60 * 60 * cycle, 14 * 24 * 60 * 60 * (cycle + 1)]
            MAlocal = solve_ivp(fun=chem, t_span=tspan, y0=self.S0,
                                args=(self.params, self.FishTank, self.B1, self.DGS), method="BDF") # if the solver does not converge try using the method "LSODA" with rtol=1e-12;atol=1e-12
            S0 = MAlocal.y[:, -1].tolist() #reintialize S0 with the last values of the previous cycle

            if cycle == 0:
                MA = MAlocal
            else:
                MA.y = np.append(MA.y, MAlocal.y, axis=1)
                MA.t = np.append(MA.t, MAlocal.t)

            n_out = Biomass(self.FishTank, MAlocal.t[-1])
            self.FishTank.Fish_number = self.FishTank.Fish_number - n_out # update the fish number
            S0[9] = self.FishTank.Fish_number * self.FishTank.Fish_weight # update the fish biomass
        # results can be store in a dataframe for further and easier exploration (e.g. TAN,TIC, CO2_removal... )
        self.df_results = pd.DataFrame(MA.y.T, columns=['CO2aq_FT', 'HCO3_FT', 'CO32_FT', 'H_FT', 'OH_FT', 'NH4_FT', 'NH3_FT',
                                                   'NO2_FT', 'Fishweight_FT', 'Fish_Biomass_FT', 'CO2aq_B1', 'HCO3_B1',
                                                   'CO32_B1', 'H_B1', 'OH_B1', 'NH4_B1', 'NH3_B1', 'NO2_B1', 'AOB_B1',
                                                   'NOB_B1', 'CO2aq_DGS', 'HCO3_DGS', 'CO32_DGS', 'H_DGS', 'OH_DGS',
                                                   'NH4_DGS', 'NH3_DGS', 'NO2_DGS'])
        self.df_results.insert(0, 'Time', MA.t)
        # all the result are in mmol/L, to convert to mg/L use the molecular weight of the species below
        self.M_NH4 = 18.04
        self.M_NH3 = 17.03
        self.M_N= 14.01
        self.M_NO2 = 46.01
        self.M_CO2 = 44.01
        self.M_HCO3= 61.01
        #save the results in a csv file
        self.df_results.to_csv('results.csv',index=False)
# alkalinity is calculated as follow Alkalinity = ([HCO3-] + 2*[CO3--] + [OH-]) * 50.04 to have it in mg/L
# TAN is calculated as follow TAN = ([NH4] + [NH3])* 14.01 to have it in mg/L
# NH3 is calculated as NH3-N in mg.l-1 the paper to be compared to experimental data : NH3-N=NH3*14.01
# Time can be divide by 60*60*24 to have it in days
# all the results in the paper are concentration in the fish tank but the results can be plot for biofilter and degasser as well

# # plot results for CO2 in the fish tank
# plt.plot(df_results['Time']/(60*60*24),df_results['CO2aq_FT']*M_CO2,label='CO2aq_FT')
# plt.xlabel('Time (days)')
# plt.ylabel('CO2 (mg/L)')
# plt.legend()
# plt.show()
#
#
# #plot alkalinity in the fish tank
# plt.plot(df_results['Time']/(60*60*24),((df_results['HCO3_FT']+2*df_results['CO32_FT']+df_results['OH_FT'])*50.04),label='Alkalinity')
# plt.xlabel('Time (days)')
# plt.ylabel('alkalinity (mg/L as CaCO3)')
# plt.legend()
# plt.show()
#
# #plot TAN in the fish tank
# plt.plot(df_results['Time']/(60*60*24),(df_results['NH4_FT']+df_results['NH3_FT'])*M_N,label='TAN')
# plt.xlabel('Time (days)')
# plt.ylabel('TAN (mg/L as N)')
#
# #twinx to plot NH3-N in mg/L
# plt.twinx()
# plt.plot(df_results['Time']/(60*60*24),df_results['NH3_FT']*M_N,label='NH3-N',color='red')
# plt.ylabel('NH3-N (mg/L)')
# plt.legend()
# plt.show()
#
# # From B1 you can plot the dosing of OH and HCO3
# B1.dosing_time = np.array(B1.dosing_time)
# B1.dosing_amount_OH = np.array(B1.dosing_amount_OH)
# B1.dosing_amount_HCO3 = np.array(B1.dosing_amount_HCO3)
# df_dosing = pd.DataFrame({'Time':B1.dosing_time/(60*60*24),'OH_dosing':B1.dosing_amount_OH,'HCO3_dosing':B1.dosing_amount_HCO3})
# df_dosing.to_csv('dosing.csv',index=False)
#
# plt.plot(B1.dosing_time/(60*60*24),B1.dosing_amount_OH,label='OH_dosing')
# plt.plot(B1.dosing_time/(60*60*24),B1.dosing_amount_HCO3,label='HCO3_dosing')
# plt.xlabel('Time (days)')
# plt.ylabel('OH and HCO3 dosing (mg/L/s)')
# plt.legend()
# plt.show()
# # pH of the FishTank
# pH = -np.log10(df_results['H_FT']*10**-3)
# plt.plot(df_results['Time']/(60*60*24),pH,label='pH')
# plt.xlabel('Time (days)')
# plt.ylabel('pH')
# plt.legend()
# plt.show()
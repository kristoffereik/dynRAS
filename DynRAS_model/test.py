from solver import Solver
from params import Params

F=60/60 #L.s-1

k1= 1.49*10**-2 #(johnson 1982) s-1
k_1 =1.89*10**4*10**-3 #(johnson 1982) L.mmol-1.s-1
k_H = 5.0*10**10*10**-3#shultz 1985 L.mmol-1.s-1
kH =5.0*10**10*10**-9.3#shultz 1985 and recalculated using K at 15C s-1
k3 =1.40*10**-3*10**3 #shultz 1985 mmol.L-1.s-1
k_3 = (1.40*10**-3/(0.45*10**-14))*10**-3 #shultz 1985 recalculated using K at 15C L.mmol-1.s-1
k_4=4.3*10**10*10**-3#eigen 1964 L.mmol-1.s-1
k4=4.3*10**10*2.73*10**-10 #eigen 1964 recalculated using K at 15C s-1

muAOB=(((0.29+0.76)/2))/(24*60*60) #s-1 pedersen 2018
muNOB=((0.28+1.04)/2)/(24*60*60)##s-1 pedersen 2018
KNH3=(1/14)*0.01#mmol.l-1 #pedersen 2018
KAlk=0.3 #mmol.l-1 #pedersen 2018
KNO2=1/14 #mmol.l-1 #pedersen 2018
rhoAOB=((0.05+0.15)/2)/(24*60*60) #s-1 pedersen 2018
rhoNOB=((0.05+0.15)/2)/(24*60*60)#s-1 pedersen 2018

YAOB=0.21*14 #gbiomass/moles-NH3 pedersen 2018
YNOB=0.04*14 #gbiomass/moles-NO2 pedersen 2018

TGC=2.7 #pedersen 2018
T=14# °C
w0=70#g

params=Params(F,k1,k_1,kH,k_H,k3,k_3,k4,k_4,muAOB,muNOB,KNH3,KAlk,KNO2,rhoAOB,rhoNOB,TGC,T,w0,YAOB,YNOB)

solver = Solver(params)
solver.solve()
import numpy as np
import pandas as pd

## Following Raymer's steps from this:
#https://github.com/jonititan/AircraftDesign/blob/master/AircraftSizing.ipynb

# Mission:
# Fly 1500nm (2778km) at M=0.6 carrying 10000lb (4535.924kg) 
# of equipment and 4 crew members at 800lb (363.874kg), 
# loiter on station for 3 hrs then return.

## Input Data:
BWB_inputs = pd.read_excel("BWB_IO.xlsx",sheet_name=1)#,usecols=["Aircraft Data","Mission Profile"])
#print(BWB_inputs)

BWB_mission = BWB_inputs.iloc[:,0:3].dropna().set_index("Mission Profile")
#print(BWB_mission)
#print(BWB_inputs["Class I"])

BWB_ClassI = BWB_inputs.iloc[:,4:7].dropna().set_index("Class I")
#print(BWB_ClassI)

BWB_Structural = BWB_inputs.iloc[:,8:11].dropna().set_index("Structural")
#print(BWB_Structural)

BWB_Aero = BWB_inputs.iloc[:,12:15].dropna().set_index("Aerodynamic")
#print(BWB_Aero)

BWB_Perf = BWB_inputs.iloc[:,16:19].dropna().set_index("Performance")
#print(BWB_Perf)

## Class I
N_pax = BWB_ClassI.loc["N_pax"][0]
N_crew = BWB_ClassI.loc["N_crew"][0]
M_pax = (N_pax + N_crew) * 80 # Average Weight of Person

## Structural:
M_cargo = BWB_Structural.loc["M_cargo"][0] # Cargo Weight [kg]
W_payload = M_pax + M_cargo

V_cr = BWB_Perf.loc["V_Cruise"][0] # Cruise Airspeed [m/s]
C_T_cruise = BWB_Perf.loc["Cruise SFC"][0] # [N/Ns] (MDO fuel consumption)
C_T_loiter = BWB_Perf.loc["Loiter SFC"][0] # [N/Ns] (MDO fuel consumption)

## Aerodynamics (Wing)
LD = BWB_Aero.loc["L/D"][0]
# segment = BWB_mission.loc["Segment"][0]
# 

## Mission Segment Fractions
mission_profile = []
for column in BWB_mission:
    segment = BWB_mission[column]
    print(segment)


"""
Range = 6000 # [km]
# Loiter Endurance [hrs -> s]
E1 = 3 * 3600 # 
E2 = 0.5 * 3600 #

# Cruise Segment:
def WF_cruise(R,C_T,V,LD):
    cruise_frac = np.exp(-( (R * C_T) / (V * LD) ))
    return cruise_frac

# Loiter Segment
def WF_loiter(E,C_T,LD):
    loiter_frac = np.exp(-( (E * C_T) / LD ))
    return loiter_frac

#Fuel Fraction
def Wf_W0(Wx_W0):
    fuel_frac = 1.06 * (1 - Wx_W0) # Extra 6% for reserve and trapped fuel
    return fuel_frac


W1_W0 = 0.97  # Warmup & TO
W2_W1 = 0.985 # Climb
W_land = 0.995 # Landing

W3_W2 = WF_cruise(Range,C_T_cruise,V_cr,LD)
W4_W3 = WF_loiter(E1,C_T_loiter,LD)
W5_W4 = W3_W2
W6_W5 = WF_loiter(E2,C_T_loiter,LD)

W7_W0 = W_land*W6_W5*W5_W4*W4_W3*W3_W2*W2_W1*W1_W0
W_FF = Wf_W0(W7_W0)

# Empty Weight Fraction Iteration (Table 3.1 &)

# def We_W0():
#     empty_frac = A * W_C0 * K_vs
#     return empty_frac

## :
A, C = 1.02, -0.06
# Variable Sweep penalty; 1.04 if A/C has variable sweep
K_vs = 1.00 # No variable sweep, thus 1.0
eps = 1.00
W0_iter = [["MTOW Guess","We/W0","We","Computed MTOW"]]
W0_i = 50e3#MTOW(WF_Empty, WF_Fuel, M_pax,M_pay)
i = 0
while i <=10:
    We_W0 = A * W0_i**C
    W0 = E1 / (1 - W_FF - We_W0)
    We = We_W0 * W0
    W0_iter.append([W0_i,We_W0,We,W0])
    eps = W0/W0_i
    W0_i = W0

    print("Computed MTOW=",W0)
    print("Eps=",eps)
    i += 1

W_fuel = W_FF * W0_i

#print(W0_iter)
print("Estimated MTOW =",W0_i,"kg")
print("Estimated Empty Weight =",We,"kg")
print("Estimated Fuel Weight =",W_fuel,"kg")
print("Payload Weight =",W_payload,"kg")


class AircraftMission:
    output_dict = {}

    def __init__(self):
        self.segments=[]
    
    def add_segment(self, seg_name, WF):
        N_seg = len(self.segments)
        self.segments.append({'segment':N_seg+1,'name':seg_name,'weight_fraction':WF})
    
    def mission_fraction(self):
        N_seg = len(self.segments)
        if N_seg <=0:
            raise ValueError("No mission segments found!")
        
        WF_mission = 1
        for i in range(N_seg):
            WF_mission = WF_mission * self.segments[i]["weight_fraction"]
        return WF_mission
    
    def W_land(self,MTOW):
        return MTOW * self.mission_fraction()
    
    def show_mission(self):
        print(self.segments)


ASW_LD_max = 16
cruise_LD_ratio = ASW_LD_max * 0.866
loiter_LD_ratio = ASW_LD_max
cruise_velocity = 596.9 # ft/s
cruise_SFC = 0.0001389
loiter_SFC = 0.0001111
mission_loiter = 10800 # 3h in s
landing_loiter = 1200 # 20 min in s
cruise_range = 9114000 # 1500nm in ft

ASW_Mission = AircraftMission()
ASW_Mission.add_segment('Warmup and Takeoff',0.97)
ASW_Mission.add_segment('Climb',0.985)


ASW_Mission.add_segment('Cruise',WF_cruise(cruise_range,cruise_SFC,cruise_velocity,cruise_LD_ratio)) 
ASW_Mission.add_segment('Loiter',WF_loiter(mission_loiter,loiter_SFC,loiter_LD_ratio)) 
ASW_Mission.add_segment('Cruise',WF_cruise(cruise_range,cruise_SFC,cruise_velocity,cruise_LD_ratio)) # 1500 nm
ASW_Mission.add_segment('Loiter',WF_loiter(landing_loiter,loiter_SFC,loiter_LD_ratio)) #20 mins


ASW_Mission.add_segment('Land',0.995)
ASW_Mission.show_mission()
"""
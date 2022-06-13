import numpy as np
import pandas as pd
from parapy.core import *

#from BWB_Inputs import BWB_Inputs

## Following Raymer's steps from this:
#https://github.com/jonititan/AircraftDesign/blob/master/AircraftSizing.ipynb

class ClassI_Estim(Base):
    #output_dict = {}

    #def __init__(self, N_pax,W_pax,W_cargo,LD_i,V_cr,SFC_cruise,SFC_loiter):
    N_pax = Input()
    N_crew = Input()
    W_pax = Input()
    W_cargo = Input()

    LD = Input()

    # LD_cr = Input(LD) * 0.866 # Cruise L/D (Raymer)
    # LD_ltr = Input(LD)          # Loiter L/D

    V_cr = Input()             # Cruise Velocity [m/s]
    SFC_cr = Input()
    SFC_ltr = Input()
    

    segments=[]
    
    @Attribute
    def segment_list(self):
        return self.segments
    
    @Attribute
    def LD_cr(self):
        return self.LD * 0.866
    
    @Attribute
    def LD_ltr(self):
        return self.LD

    def add_segment(self, seg_name, value,info=None):
        N_seg = len(self.segments)

        if seg_name=="Cruise":
            WF_cr = self.WF_cruise(value,self.SFC_cr,self.V_cr,self.LD_cr)
            self.segments.append({'segment':N_seg+1,'name':seg_name,'weight_fraction':WF_cr,'info':info})
        
        elif seg_name == "Loiter":
            WF_ltr = self.WF_loiter(value,self.SFC_cr,self.LD_cr)
            self.segments.append({'segment':N_seg+1,'name':seg_name,'weight_fraction':WF_ltr,'info':info})
        
        else: # Assume fraction was given
            self.segments.append({'segment':N_seg+1,'name':seg_name,'weight_fraction':value,'info':info})

    # Cruise Segment Fraction
    def WF_cruise(self,R,C_T,V,LD):
        cruise_frac = np.exp(-( (R*1e3 * C_T) / (V * LD) ))
        return cruise_frac

    # Loiter Segment Fraction
    def WF_loiter(self,E,C_T,LD):
        loiter_frac = np.exp(-( (E*3600 * C_T) / LD ))
        return loiter_frac
    
    @Attribute
    def W_F_loiter(self):
        return self.WF_loiter
    
    @Attribute
    def W_F_cruise(self):
        return self.WF_cruise

    # Total Mission Fraction
    @Attribute
    def Wx_W0(self):
        N_seg = len(self.segments)
        if N_seg <=0:
            raise ValueError("No mission segments found!")
        
        WF_mission = 1
        for i in range(N_seg):
            WF_mission = WF_mission * self.segments[i]["weight_fraction"]
        print("Wx_W0=",WF_mission)
        return WF_mission
    
    @Attribute
    def W_FF(self): # Fuel Fraction
        # Extra 6% for reserve and trapped fuel (Raymer)
        fuel_frac = 1.06 * (1 - self.Wx_W0) 
        return fuel_frac

    @Attribute
    def MTOW(self):
        # MTOW eq according to Raymer
        A, C = 1.02, -0.06
        # Variable Sweep penalty; 1.04 if A/C has variable sweep
        K_vs = 1.00 # No variable sweep, thus 1.0
        
        W0_iter = [["MTOW Guess","We/W0","We","Computed MTOW"]]
        W0_i = 600e3 # Initial Guess
        i = 0
        eps = 1.00
        while i <=20:
            We_W0 = A * (W0_i**C)
            W0 = (self.W_cargo + (self.W_pax*(self.N_pax+self.N_crew))) / (1 - self.W_FF - We_W0)
            print("W0 =",W0)
            print("We_W0 =",We_W0)
            We = We_W0 * W0
            W0_iter.append([W0_i,We_W0,We,W0])
            eps = W0/W0_i
            W0_i = W0

            #print("Computed MTOW=",W0)
            #print("Eps=",eps)
            i += 1

        return W0_i, We
    @Attribute
    def W_end(self):
        return self.MTOW[0] * self.Wx_W0
    
    @Attribute
    def W_fuel(self):
        return self.W_FF * self.MTOW[0]

    # @Attribute
    # def show_mission(self):
    #     # print("Estimated MTOW =",self.MTOW,"kg")
    #     # print("Estimated Empty Weight =",,"kg")
    #     # print("Estimated Fuel Weight =",W_fuel,"kg")
    #     # print("Payload Weight =",W_payload,"kg")
    #     print(self.segments)

################ Retrieve Input Data from XLSX: ##########################

# ## Input Data:
# BWB_inputs = pd.read_excel("BWB_IO.xlsx",sheet_name=1)#,usecols=["Aircraft Data","Mission Profile"])
# #print(BWB_inputs)

# BWB_mission = BWB_inputs.iloc[:,0:3].dropna().set_index("Mission Profile")
# #print(BWB_mission)
# #print(BWB_inputs["Class I"])

# BWB_ClassI = BWB_inputs.iloc[:,4:7].dropna().set_index("Class I")
# #print(BWB_ClassI)

# BWB_Structural = BWB_inputs.iloc[:,8:11].dropna().set_index("Structural")
# #print(BWB_Structural)

# BWB_Aero = BWB_inputs.iloc[:,12:15].dropna().set_index("Aerodynamic")
# #print(BWB_Aero)

# BWB_Perf = BWB_inputs.iloc[:,16:19].dropna().set_index("Performance")
# #print(BWB_Perf)

# ## Class I
# N_pax = BWB_ClassI.loc["N_pax"][0]
# N_crew = BWB_ClassI.loc["N_crew"][0]
# M_pax = (N_pax + N_crew) * 80 # Average Weight of Person

# ## Structural:
# M_cargo = BWB_Structural.loc["M_cargo"][0] # Cargo Weight [kg]
# W_payload = M_pax + M_cargo

# V_cr = BWB_Perf.loc["V_Cruise"][0] # Cruise Airspeed [m/s]
# C_T_cruise = BWB_Perf.loc["Cruise SFC"][0] # [N/Ns] (MDO fuel consumption)
# C_T_loiter = BWB_Perf.loc["Loiter SFC"][0] # [N/Ns] (MDO fuel consumption)

# ## Aerodynamics (Wing)
# LD_i = BWB_Aero.loc["L/D"][0]

# ## Mission Reader:

# ## Test AC Mission Class
# BWB_Mission = AircraftMission(N_pax,M_pax,M_cargo,LD_i,V_cr,C_T_cruise,C_T_loiter)

# for index in range(len(BWB_mission)):
    
#     segment = BWB_mission.index.values[index]
#     value = BWB_mission.iloc[index]["Value"]
#     unit = BWB_mission.iloc[index]["Unit"]
#     #print([segment,value,unit])
#     BWB_Mission.add_segment(segment,value,unit)

# BWB_Mission.show_mission()
# print("Fuel Weight=",BWB_Mission.W_fuel())
# BWB_Mission.MTOW()
# BWB_Mission.W_end()


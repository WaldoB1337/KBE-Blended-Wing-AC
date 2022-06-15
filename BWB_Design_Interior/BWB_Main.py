from turtle import position
from parapy.core import *
from parapy.geom import *

import pandas as pd
import os

from BWB_Class_I import ClassI_Estim
from BWB_CentralFuselage import CentralFuselage
from parapy.exchange.step import STEPWriter
################ Retrieve Input Data from XLSX: ##########################
DIR = os.path.dirname(__file__)
## Input Data:
BWB_inputs = pd.read_excel("BWB_IO.xlsx",sheet_name=1)#,usecols=["Aircraft Data","Mission Profile"])
#print(BWB_inputs)

BWB_segments = BWB_inputs.iloc[:,0:3].dropna().set_index("Mission Profile")
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
M_pax = 80 # [kg] Average Weight of Person

## Structural:
M_cargo = BWB_Structural.loc["M_cargo"][0] # Cargo Weight [kg]
W_payload = M_pax + M_cargo

V_cr = BWB_Perf.loc["V_Cruise"][0] # Cruise Airspeed [m/s]
C_T_cruise = BWB_Perf.loc["Cruise SFC"][0] # [N/Ns] (MDO fuel consumption)
C_T_loiter = BWB_Perf.loc["Loiter SFC"][0] # [N/Ns] (MDO fuel consumption)

## Aerodynamics (Wing)
LD_i = BWB_Aero.loc["L/D"][0]

class BWB_Design(GeomBase):
    N_pax = Input()
    N_crew = Input()
    #self.pax_dist = pax_dist
    LD = Input()
    V_cr = Input()      # Cruise Velocity [m/s]
    SFC_cruise = Input()
    SFC_loiter = Input()
    W_pax = Input()
    W_cargo = Input()
    #W_fuel = Input()

    @Attribute
    def Mission(self):
        BWB_mission= ClassI_Estim(N_pax=self.N_pax,N_crew=self.N_crew,
                                W_pax=M_pax,W_cargo=self.W_cargo,
        LD=self.LD,V_cr=self.V_cr,SFC_cr=self.SFC_cruise,SFC_ltr=self.SFC_loiter)

        for index in range(len(BWB_segments)):
            segment = BWB_segments.index.values[index]
            value = BWB_segments.iloc[index]["Value"]
            unit = BWB_segments.iloc[index]["Unit"]
            #print([segment,value,unit])
            BWB_mission.add_segment(segment,value,unit)
        
        return BWB_mission

    @Part(parse=False)
    def Interior(self):
        #print(self.Mission.W_fuel)
        interior = CentralFuselage(N_pax=self.Mission.N_pax,N_crew=self.Mission.N_crew, 
                            m_cargo=self.Mission.W_cargo,W_fuel=self.Mission.W_fuel)
        # interior.orientation.rotate("z",-90,deg=True)
        # interior.position.rotate("z",-90,deg=True)
        return interior
        # return CentralFuselage(N_pax=self.N_pax,N_crew=self.N_crew, 
        #                     m_cargo=self.W_cargo,W_fuel=self.Mission.W_fuel)
    
    # @Part
    # def Wing(self):
    #     return None
    
    # @Part
    # def Engines(self):
    #     return None

    @Attribute
    def CAD_Export(self):
        print(DIR)
        return STEPWriter(trees=[self.Interior.cabin, self.Interior.cargo,
                                self.Interior.fuel_tanks],
                          default_directory=DIR,
                          filename="BWB_CAD.stp")
    

if __name__ == '__main__':
    from parapy.gui import display

    obj = BWB_Design(N_pax=N_pax,N_crew=N_crew,W_pax=M_pax,
                    W_cargo=M_cargo,LD=LD_i,V_cr=V_cr,
                    SFC_cruise=C_T_cruise,SFC_loiter=C_T_loiter)
    display(obj)
    #display(obj.ellipses)
    obj.CAD_Export.write()
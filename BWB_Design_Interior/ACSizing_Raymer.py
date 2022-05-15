#https://github.com/jonititan/AircraftDesign/blob/master/AircraftSizing.ipynb

# Mission:
# Fly 1500nm (2778km) at M=0.6 carrying 10000lb (4535.924kg) 
# of equipment and 4 crew members at 800lb (363.874kg), 
# loiter on station for 3 hrs then return.

import matplotlib
import matplotlib.pyplot as plt
import numpy

## Take-Off Estimation
def TGW(EmptyW_F, Fuel_F, Crew_W, Payload_W):
    return (Crew_W+Payload_W)/(1-Fuel_F-EmptyW_F)
#e.g. ASW example
Crew_Weight = 800 #lb
Payload_Weight = 10000 #lb
Fuel_Fraction = 0.387
Empty_Weight_Fraction = 0.4361
print("Crew Weight %.2f lb" % (Crew_Weight))
print("Payload Weight %.2f lb" % (Payload_Weight))
print("Fuel Fraction %.4f" % (Fuel_Fraction,))
print("Empty Weight Fractiont %.4f" % (Empty_Weight_Fraction,))
TGW_tmp_1 = TGW(Empty_Weight_Fraction, Fuel_Fraction, Crew_Weight, Payload_Weight)
print("Takeoff Gross Weight W_0 %.2f lb" % (TGW_tmp_1))

## Empty Weight Estimation:
def EWF(TakeoffG_W,HS_metric,WeightS_E,VWS=False,VWS_penalty=1.04):
    if VWS:
        VWS_p_factor = VWS_penalty
    else:
        VWS_p_factor = 1
        
    return HS_metric*pow(TakeoffG_W, WeightS_E)*VWS_p_factor

bomber_HS_metric = 0.93
bomber_weight_sensitivty_e = -0.07
estimated_TGW = 50000
print(EWF(estimated_TGW,bomber_HS_metric,bomber_weight_sensitivty_e))

## Fuel Fraction 
def FF(mission_ratio,RT_Allowance=1.06):
    return RT_Allowance*(1-mission_ratio)

## Cruise Fuel 
def cruise_WF(Range, Specific_FC, velocity, LD_ratio):
    return numpy.exp((-1*Range*Specific_FC)/(velocity*LD_ratio))

high_bypass_turbofan_cruise_SFC = 0.0001389

## Loiter Fuel Burn

def loiter_WF(Endurance, Specific_FC, LD_ratio):
    return numpy.exp((-1*Endurance*Specific_FC)/(LD_ratio))

high_bypass_turbofan_loiter_SFC = 0.0001111

##L/D Estimation:
class AircraftMission:
    variable_stan_disp = {'name':{'description':'Name:','format':'%s'},
                         'weight_fraction':{'description':'WF:','format':'%.4f'},
                         'description':{'description':'Description:','format':'%s'},
                         'segment':{'description':'Segment:','format':'%i'}}
    def __init__(self):
        self.segments=[]
    
    def append_segment(self, seg_name, seg_WF, seg_description='None Provided'):
        no_of_segments = len(self.segments)
        self.segments.append({'segment':no_of_segments+1,'name':seg_name,'weight_fraction':seg_WF,'description':seg_description})
    
    def mission_weight_ratio(self):
        no_of_segments = len(self.segments)
        if no_of_segments <= 0:
            raise ValueError('No mission segments have been added to the mission!')
        
        mwr = 1
        for a in range(0,no_of_segments):
            mwr = mwr*self.segments[a]['weight_fraction']
        return mwr
    
    def landing_weight(self, TakeoffG_W):
        return TakeoffG_W * self.mission_weight_ratio()
    
    def display_mission(self,itp=['segment','name','weight_fraction','description'],p_dict=variable_stan_disp):
        no_of_segments = len(self.segments)
        if no_of_segments <= 0:
            raise ValueError('No mission segments have been added to the mission!')
        
        for a in range(0,no_of_segments):
            disp_list = [p_dict[item_a]['description']+p_dict[item_a]['format']+' ' for item_a in itp]
            variable_list = [self.segments[a][item_b] for item_b in itp]
            disp_string = ''.join(disp_list) % tuple(variable_list)
            print(disp_string)

ASW_LD_max = 16
cruise_LD_ratio = ASW_LD_max * 0.866
loiter_LD_ratio = ASW_LD_max
cruise_velocity = 596.9 # ft/s
cruise_SFC = high_bypass_turbofan_cruise_SFC
loiter_SFC = high_bypass_turbofan_loiter_SFC
mission_loiter = 10800 # 3h in s
landing_loiter = 1200 # 20 min in s
cruise_range = 9114000 # 1500nm in ft

ASW_Mission = AircraftMission()
ASW_Mission.append_segment('Warmup and Takeoff',0.97)
ASW_Mission.append_segment('Climb',0.985)
ASW_Mission.append_segment('Cruise',cruise_WF(cruise_range,cruise_SFC,cruise_velocity,cruise_LD_ratio),seg_description='1500 nm range') 
ASW_Mission.append_segment('Loiter',loiter_WF(mission_loiter,loiter_SFC,loiter_LD_ratio),seg_description='3 hr Loiter') 
ASW_Mission.append_segment('Cruise',cruise_WF(cruise_range,cruise_SFC,cruise_velocity,cruise_LD_ratio),seg_description='1500 nm range') # 1500 nm
ASW_Mission.append_segment('Loiter',loiter_WF(landing_loiter,loiter_SFC,loiter_LD_ratio),seg_description='20 minute loiter') #20 mins
ASW_Mission.append_segment('Land',0.995)
ASW_Mission.display_mission()

MWR_tmp = ASW_Mission.mission_weight_ratio()
print("Misson Weight Ratio %.4f" % (MWR_tmp,))
FF_tmp = FF(MWR_tmp)
print("Fuel Fraction %.4f" % (FF_tmp,))
TakeoffG_W_estimate = 50000
EWF_tmp = EWF(TakeoffG_W_estimate,0.93,-0.07)
print("Empty Weight Fraction %.5f based on TGW W_0 estimate %.3f lb" % (EWF_tmp,TakeoffG_W_estimate))
TGW_tmp = TGW(EWF_tmp , FF_tmp, Crew_Weight, Payload_Weight)
print("Takeoff Gross Weight W_0 %.3f lb"  % (TGW_tmp))
LGW_tmp = ASW_Mission.landing_weight(TGW_tmp)
print("Landing Gross Weight W_x %.3f lb"  % (LGW_tmp))
fuel_burn = TGW_tmp-LGW_tmp
print("Fuel Burned %.3f lb"  % (fuel_burn))


    
def calc_WG0(WGV,Fuel_F,Crew_W, Payload_W, HS_metric,WeightS_E,VWS=False,VWS_penalty=1.04):
    EWF_calc = EWF(WGV,HS_metric,WeightS_E,VWS=False,VWS_penalty=1.04)
    return TGW(EWF_calc, Fuel_F, Crew_W, Payload_W)

def diff_calc_WG0(WGV,Fuel_F,Crew_W, Payload_W, HS_metric,WeightS_E,VWS=False,VWS_penalty=1.04):
    tmp_wgv = calc_WG0(WGV,Fuel_F,Crew_W, Payload_W, HS_metric,WeightS_E,VWS,VWS_penalty)
    result = WGV - tmp_wgv
    if result < 0:
        mag_result = result * -1
    else:
        mag_result=result
    return mag_result
    

def calc_WG0_array(WGV_list,Fuel_F,Crew_W, Payload_W, HS_metric,WeightS_E,VWS=False,VWS_penalty=1.04):
    WG0_array=[]
    EWF_array=[]
    diff_array=[]
    for a in range(0,len(WGV_list)):
        WG0_array.append(calc_WG0(WGV_list[a],Fuel_F,Crew_W, Payload_W, HS_metric,WeightS_E,VWS,VWS_penalty))
        EWF_array.append(EWF(WGV_list[a],HS_metric,WeightS_E,VWS=False,VWS_penalty=1.04))
        diff_array.append(diff_calc_WG0(WGV_list[a],Fuel_F,Crew_W, Payload_W, HS_metric,WeightS_E,VWS,VWS_penalty))
    return WG0_array, EWF_array, diff_array

WGV_array = [50000]
for a in range(0,250):
    WGV_array.append(WGV_array[a]+50)
    
calc_WG0_tmp, calc_EWF_tmp, diff_array = calc_WG0_array(WGV_array,FF_tmp,Crew_Weight, Payload_Weight, bomber_HS_metric,bomber_weight_sensitivty_e)
best_guess_index = numpy.argmin(diff_array)
fig, ax = plt.subplots()
axes = [ax, ax.twinx(), ax.twinx()]
fig.subplots_adjust(right=0.75)
axes[-1].spines['right'].set_position(('axes', 1.2))
axes[0].plot(WGV_array,calc_WG0_tmp,c='b')
axes[0].set_title('Estimated vs Calculated Takeoff Gross Weight for\nASW Mission Best Guess $W_0$:%.3f lb EWF:%.4f' % (WGV_array[best_guess_index],calc_EWF_tmp[best_guess_index]))
axes[0].set_ylabel('Calculated $W_0$ (lb)',color='b')
axes[0].set_xlabel('Estimated $W_0$ (lb)')
axes[0].tick_params('y', colors='b')
axes[1].plot(WGV_array,calc_EWF_tmp,c='r')
axes[1].set_ylabel('Empty Weight Fraction',color='r')
axes[1].tick_params('y', colors='r')
axes[2].plot(WGV_array,diff_array,c='g')
axes[2].set_ylabel('|difference|',color='g')
axes[2].tick_params('y', colors='g')
plt.show()

## optimization:


import scipy.optimize as optimize
minimum_W_0 = optimize.minimize(diff_calc_WG0,50000, args = (FF_tmp,Crew_Weight, Payload_Weight, bomber_HS_metric,bomber_weight_sensitivty_e))
min_EWF = EWF(minimum_W_0.x,0.93,-0.07)
print("Optimised Takeoff Gross Weight W_0 found %.3f lb with EWF %.4f" % (minimum_W_0.x,min_EWF))


print('Difference between the Takeoff Gross Weight W_0 estimations: {} lbs'.format((WGV_array[best_guess_index] - minimum_W_0.x)[0]))
print('Difference between the Empty Weight Fraction EWF estimations: {} %'.format(100*((calc_EWF_tmp[best_guess_index] - min_EWF)[0])))
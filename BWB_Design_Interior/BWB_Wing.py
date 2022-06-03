from parapy.core import *
from parapy.geom import *
import numpy as np
from kbeutils.geom.curve import Naca5AirfoilCurve, Naca4AirfoilCurve
import kbeutils.avl as avl
#import parapy.lib.avl as avl

from BWB_Cabin import Cabin
from BWB_Fuel import FuelTank


class Wing(GeomBase):
    taper = Input(0.55)
    surface_scale=Input(65)
    airfoil_name = Input(str(2312))
    cabin = Input()
    fuel_tanks = Input()
    origin = Input()
    wingspan = Input()

    # def __init__(self,coor,cabin,fuel_tanks):
    #     self.origin = coor
    #     self.cabin = cabin
    #     self.fuel_tanks = fuel_tanks
        

    @Attribute#(parse=False)
    def AirfoilCenter(self):
        Airfoil = DynamicType(type=Naca5AirfoilCurve \
                           if len(self.airfoil_name) == 5 else Naca4AirfoilCurve,
                           designation=self.airfoil_name,
                           mesh_deflection=0.00001,
                           hidden=True)

        airfoil_curve = ScaledCurve(curve_in=Airfoil,
                           reference_point=XOY,
                           factor=self.surface_scale)#2.5*self.cabin.l_cabin_tot)
 
        return airfoil_curve
        # TransformedCurve(airfoil_curve,
        #                     from_position=XOY.rotate("z",-90,deg=True),
        #                     to_position=self.origin_front[1].position.translate(
        #                                             "y",-0.1*self.cabin.l_cabin_tot))

    @Attribute
    def CabinShells(self):
        """Potentiall add more shells along a line for resolution or Bezier ?"""
        origin = self.origin_front[1].position.translate("y",-0.1*self.cabin.l_cabin_tot)
        #origin = XOY.translate("y",-0.1*self.cabin.l_cabin_tot)
        print(origin)
        airfoil_shells = []
        tank_faces = [self.fuel_tanks.tank_tip1,
                    self.fuel_tanks.tank_tip2]
        
        airfoil_points = [origin]

        # Cabin Airfoils at front control points
        for location in self.cabin.control_points:
            #print(location)
            point = location
            airfoil_points.append([point.x,point.y,point.z])
            
        # Fuel Tank Airfoils at Chord Face locations
        for face in tank_faces:
            point = face.position.location
            airfoil_points.append([point.x,point.y,point.z])

        airfoil_points_s = sorted(airfoil_points,key=lambda l:l[0])
        
        for i in range(len(airfoil_points_s)):
            loc = airfoil_points_s[i]
            point = Point(loc[0],loc[1],loc[2])
            #print(point)
            if i==0 or i==len(airfoil_points_s)-1:
                shell_i= ScaledCurve(curve_in=self.AirfoilCenter,
                           reference_point=XOY,
                           factor=self.taper) # Taper Ratio + 5% Margin ?
            
                shell_i = TransformedCurve(shell_i,
                                from_position=XOY.rotate("z",-90,deg=True),
                                to_position=point)
                shift = 0.8 * (shell_i.bbox.center -shell_i.position)
                #print(shift)
                shell = TranslatedCurve(curve_in=shell_i,
                                        displacement=-shift)
            else: 
                shell = TransformedCurve(self.AirfoilCenter,
                            from_position=XOY.rotate("z",-90,deg=True),
                            to_position=point)
            airfoil_shells.append(shell)

        #print(airfoil_points_s)
        return airfoil_shells, airfoil_points_s

    @Part(parse=False)
    def WingShells(self):
        #airfoil = 5
        reso = 10
        taper_init = 0.75* self.fuel_tanks.tank_c_r / (0.75*self.cabin.l_cabin)
        print(taper_init)
        taper_line = np.linspace(taper_init,self.taper,reso)
        
        tank_airfoils = []
        ref_points = [self.fuel_tanks.fuel_tank1.profile1.center,
                    self.fuel_tanks.fuel_tank2.profile1.center]
        for ref in range(len(ref_points)):
            line_vec = self.fuel_tanks.midline[ref] / reso
            
            for i in range(reso):
                point = ref_points[ref] + (i+1)*line_vec
                #print(point)
                shell_i= ScaledCurve(curve_in=self.AirfoilCenter,
                           reference_point=XOY,
                           factor=taper_line[i]) # Taper Ratio + 5% Margin ?
            
                shell_i = TransformedCurve(shell_i,
                                from_position=XOY.rotate("z",-90,deg=True),
                                to_position=point)

                shift = 0.8 * (shell_i.bbox.center -shell_i.position)
                #print(shift)
                shell = TranslatedCurve(curve_in=shell_i,
                                        displacement=-shift)
                tank_airfoils.append(shell)
        return tank_airfoils

    @Attribute
    def WingletShells(self):
        height = 5
        winglet_angle = 45
        winglet_shells = []

        frame_1 = Point(self.CabinShells[1][-1][0],
                        self.CabinShells[1][-1][1],
                        self.CabinShells[1][-1][2])\
                            .translate("z",height)#.rotate("z",-90,deg=True)
        frame_2 = Point(self.CabinShells[1][0][0],
                        self.CabinShells[1][0][1],
                        self.CabinShells[1][0][2])\
                            .translate("z",height)#.rotate("z",-90,deg=True)
        ref_point = [frame_1,frame_2]
        for ref in range(len(ref_point)):
            if ref==0:
                k = 1
            elif ref==1:
                k = -1
            point = ref_point[ref] + Vector(k*height*np.cos((np.pi/180)*winglet_angle),
                                            0,height*np.sin((np.pi/180)*winglet_angle))
            
            shell= ScaledCurve(curve_in=self.AirfoilCenter,
                            reference_point=XOY,
                            factor=self.taper*0.5) # Taper Ratio + 5% Margin ?
                
            shell = TransformedCurve(shell,
                            from_position=XOY.rotate("z",-90,deg=True),
                            to_position=point)
            
            shift = 0.8 * (shell.bbox.center -shell.position)
            shell = TranslatedCurve(curve_in=shell,
                                        displacement=-shift)

            shell = RotatedCurve(curve_in=shell,
                            rotation_point=shell.position,
                            vector=Vector(0, 1, 0),
                            angle=np.radians(-k*45))
            winglet_shells.append(shell)

        return winglet_shells
    
    @Part(parse=False)
    def WingletFrame(self):
        return self.WingletShells
    
    @Part(parse=False)
    def CabinFrame(self):
        #""" ADD WINGLETS AT ENDS OF THE LIST"""
        # shells = self.WingShells[0]
        # # Left Side:
        # shells.insert(0,self.WingletShells[1])
        # # Right Side:
        # shells.append(self.WingletShells[0])
        return self.CabinShells[0]
    
    @Attribute
    def CombinedFrame(self):
        shells = self.CabinFrame
        # Left Side:
        shells.insert(0,self.WingletShells[1])
        # Right Side:
        shells.append(self.WingletShells[0])
        return shells
    
    @Part(parse=False)
    def LiftSurface(self):
        # shells = self.CabinFrame
        # # Left Side:
        # shells.insert(0,self.WingletShells[1])
        # # Right Side:
        # shells.append(self.WingletShells[0])
        return RuledSolid(profiles=self.CombinedFrame)
    
    @Attribute
    def wing_area(self):
        return 0.5 *self.surface_scale*(1+self.taper) * self.wingspan
    
    @Attribute
    def wing_span(self):
        return self.wingspan
    
    @Attribute
    def wing_mac(self):
        return 2/3*self.surface_scale*(1+self.taper+self.taper**2)/(1+self.taper)
    
    @Attribute
    def wing_cog(self):
        return self.LiftSurface.cog
    
    @Attribute
    def root_section(self):
        # print("Root Section")
        # print(self.CabinFrame[3])
        return self.CabinFrame[3]

    @Attribute
    def tip_section(self):
        # print("Tip Section")
        # print(self.WingShells[9])
        return self.WingShells[9]
    
    @Part(parse=False)
    def avl_root(self):
        # mid = int(0.5*len(self.CombinedFrame))
        # half_section = self.CombinedFrame[:mid]
        return avl.SectionFromCurve(curve_in=self.root_section)
    
    @Part(parse=False)
    def avl_tip(self):
        # mid = int(0.5*len(self.CombinedFrame))
        # half_section = self.CombinedFrame[:mid]
        return avl.SectionFromCurve(curve_in=self.tip_section)
    
    @Part
    def avl_Surface(self):
        return avl.Surface(name="BWB_Lifting_Surface",
                           n_chordwise=self.surface_scale,
                           chord_spacing=avl.Spacing.cosine,
                           n_spanwise=self.wing_span,
                           span_spacing=avl.Spacing.cosine,
                           y_duplicate= 0,#self.position.point[1],#if self.is_mirrored else None,
                           sections=[self.avl_root,self.avl_tip])

if __name__ == '__main__':
    from parapy.gui import display

    obj = Wing(XOY,Cabin,FuelTank)
    display(obj)
    #display(obj.ellipses)
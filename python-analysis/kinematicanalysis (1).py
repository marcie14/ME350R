# -*- coding: utf-8 -*-
"""
Created on Sun 31 Mar

@author: Meredith, Alexa
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numerical

class FourBarVisualize():
    
    # Constructor
    def __init__(self, length = 2.0):
        self.l_vec = [2.23, 3, 3.94, 3]    # length of link 2
        self.omega_2 = -1.0     # constant rad/sec clockwise
        self.offset = np.rad2deg(np.arctan2(0.4476,1.821))
        
        self.fig, self.ax = plt.subplots(1,1,figsize=(6,6))   # New figure
        self.ax.set_aspect('equal')
        self.fig.set_tight_layout(True)
        self.ax.set_xlim([-3*length,3*length])
        self.ax.set_ylim([-3*length,3*length])
    
    def fourBarPlot(self,th_2_in,showVelArrows=True):
        self.ax.clear()
        
        th_pos_soln = self.fourbarpos(0,th_2_in,-1)
        
        th_pos_soln = [self.getBoundedAngle(x+self.offset) for x in th_pos_soln]
        (omega_vec_out,VA,VBA,VB) = self.fourbarvel(th_pos_soln,self.omega_2)
        
        # Get position of A with respect to origin at O2 (GCS)
        th_pos_soln = [np.deg2rad(x) for x in th_pos_soln]
        R_AO2X = self.l_vec[1]*np.cos(th_pos_soln[1])
        R_AO2Y = self.l_vec[1]*np.sin(th_pos_soln[1])

        # Get position of B with respect to A
        R_BAX = self.l_vec[2]*np.cos(th_pos_soln[2])
        R_BAY = self.l_vec[2]*np.sin(th_pos_soln[2])

        # Get position of B in the GCS
        R_BX = R_AO2X + R_BAX
        R_BY = R_AO2Y + R_BAY

        # Get position of O4 in GCS
        R_O4X = self.l_vec[0]*np.cos(th_pos_soln[0])
        R_O4Y = self.l_vec[0]*np.sin(th_pos_soln[0])
        
        self.ax.set_xlim([-8,8])
        self.ax.set_ylim([-8,8])
        self.ax.set_aspect('equal')
        
        self.ax.plot([0,R_AO2X],[0,R_AO2Y],color='red', linestyle='solid')                       # Link 2
        self.ax.scatter([R_AO2X],[R_AO2Y], s=200, marker='o', facecolors='none', edgecolors='r') # Joint A
        self.ax.scatter([0],[0], s=200, marker='o', facecolors='none', edgecolors='r')     # Joint O_2
        
        self.ax.plot([R_AO2X,R_BX],[R_AO2Y,R_BY],color='blue', linestyle='solid')                       # Link 3
        self.ax.scatter([R_BX],[R_BY], s=200, marker='o', facecolors='none', edgecolors='b') # Joint B
        
        self.ax.plot([R_O4X,R_BX],[R_O4Y,R_BY],color='green', linestyle='solid')                       # Link 4
        self.ax.scatter([R_O4X],[R_O4Y], s=200, marker='o', facecolors='none', edgecolors='g') # Joint O4
        self.ax.plot([0,R_O4X],[0,R_O4Y],color='black', linestyle='solid')                       # Link 1

        # If we want to see the velocity vector (with arrows!)
        if showVelArrows:
            self.ax.arrow(R_AO2X,R_AO2Y,VA[0],VA[1],head_width=0.1,width=0.001,length_includes_head=True,linestyle='dashed',color='red',overhang=1.0)
            self.ax.arrow(R_BX,R_BY,VBA[0],VBA[1],head_width=0.1,width=0.001,length_includes_head=True,linestyle='dashed',color='blue',overhang=1.0)
            self.ax.arrow(R_BX,R_BY,VB[0],VB[1],head_width=0.1,width=0.001,length_includes_head=True,linestyle='dashed',color='green',overhang=1.0)
   
        
        # Centers the axes
        self.ax.spines['left'].set_position('center') 
        self.ax.spines['bottom'].set_position('center')
        self.ax.spines['right'].set_color('none')
        self.ax.spines['top'].set_color('none')

    def fourbarpos(self,th_1,th_2,delta):
        
        # K constants based on link lengths
        K1 = self.l_vec[0]/self.l_vec[1]
        K2 = self.l_vec[0]/self.l_vec[3]
        K3 = (self.l_vec[1]**2 - self.l_vec[2]**2 + self.l_vec[3]**2 + self.l_vec[0]**2)/(2 * self.l_vec[1] * self.l_vec[3])
        K4 = self.l_vec[0]/self.l_vec[2]
        K5 = (self.l_vec[3]**2 - self.l_vec[0]**2 - self.l_vec[1]**2 - self.l_vec[2]**2)/(2 * self.l_vec[1] * self.l_vec[2])
        
        
        # Quadratic constants
        th_1 = np.deg2rad(th_1)
        th_2 = np.deg2rad(th_2)  # converting to rads
        A = -K1 - K2*np.cos(th_2) + K3 + np.cos(th_2)
        B = -2*np.sin(th_2)
        C = K1 - K2*np.cos(th_2) + K3 - np.cos(th_2)
        D = np.cos(th_2) - K1 +K4*np.cos(th_2) + K5
        E = -2*np.sin(th_2)
        F = K1 + (K4 - 1)*np.cos(th_2) + K5
        
        # We don't want to handle imaginary roots
        disc_4 = B**2-(4*A*C)
        disc_3 = E**2-(4*D*F)
        if disc_4 < 0 or disc_3 < 0:
            print (B,A,C,E,D,F)
            print (disc_3,disc_4)
            print('rip2')
            raise SystemExit('Error: This function does not handle imaginary roots')
        
        # Solve for thetas 
        th_4 = 2*np.arctan2((-B + delta*np.sqrt(B**2-4*A*C)),(2*A))
        th_3 = 2*np.arctan2((-E + delta*np.sqrt(E**2-4*D*F)),(2*D))

        #print(th_1,th_2,th_3,th_4)

        return [np.rad2deg(x) for x in [th_1,th_2,th_3,th_4]]
    
    def fourbarvel(self,th_vec,omega_2):
        th_vec = [np.deg2rad(x) for x in th_vec]
        #print(np.rad2deg(th_vec[2]),th_vec[3],th_vec[2]-th_vec[3])
        # solve for the angular velocities
        omega_3 = omega_2*self.l_vec[1]/self.l_vec[2]*(np.sin(th_vec[3]-th_vec[1]))/(np.sin(th_vec[2]-th_vec[3]))
        try:
            omega_4 = omega_2*self.l_vec[1]/self.l_vec[3]*(np.sin(th_vec[1]-th_vec[2]))/(np.sin(th_vec[3]-th_vec[2]))
        except ZeroDivisionError:
            print('rip')
            omega_4 = float(0)
        except Exception as e:
            print(f"An error occurred: {e}")
            omega_4 = float(0)
        # solve for the absolute velocities
        VA = [self.l_vec[1]*omega_2*x for x in [-np.sin(th_vec[1]),np.cos(th_vec[1])]]
        VBA = [self.l_vec[2]*omega_3*x for x in [-np.sin(th_vec[2]),np.cos(th_vec[2])]]
        VB = [self.l_vec[3]*omega_4*x for x in [-np.sin(th_vec[3]),np.cos(th_vec[3])]]

        omega_vec_out = [0,omega_2,omega_3,omega_4]
        return (omega_vec_out,VA,VBA,VB)
        
    def getBoundedAngle(self,angle):
        if angle < 0:
            return 360+angle
        elif angle > 360:
            return angle-360
        else:
            return angle
    
    # Function to animate the crank
    def animate(self):
             
        # Runs on the first animation frame
        def init():
            self.ax.clear()
            
        # Adds the four bar plotting function to the animation 
        def getFrames(frames):
            self.fourBarPlot(frames-self.offset)
            
        
        # Animation function
        #delay_interval = round(1/abs(self.omega_2*180/np.pi)*2*1000)
        delay_interval = 35
        theta_frames = np.linspace(0,361,361)
        anim = FuncAnimation(self.fig, getFrames, init_func=init, frames=theta_frames[::-1], interval=delay_interval,repeat_delay=100)
    
        # Save as GIF
        anim.save('anim2.gif', dpi=200, writer='pillow') 
    
    
if __name__ == '__main__':
    l_2 = 2.0
    #mechanism = FourBarVisualize(l_2) # Creates an instance with link 2 of length 2.0
    hy = np.hypot(1.821, 0.4476)
    mechanism = FourBarVisualize(l_2)
    mechanism.l_vec = [1.875, 3, 4, 3.06]
    th1, th2, th3, th4 = mechanism.fourbarpos(0,np.linspace(0,361,1),-1)
    (omega_vec_out,VA,VBA,VB) = mechanism.fourbarvel([th1, th2, th3, th4], -1)
    mechanism.offset = np.rad2deg(np.arctan2(0.4476,1.821))
    mechanism.animate()

    
    mechanism2 = FourBarVisualize()
    mechanism2.l_vec = [hy,3.07181611,3.99108686,3.05088475]
    print('b')
    #print('th4',th1, th2, th3, th4)
    th1_2, th2_2, th3_2, th4_2 = mechanism2.fourbarpos(0,th4 - 39.01,-1)
    #print(th1, th2, th3, th4)
    mechanism2.offset = - np.rad2deg(np.arctan2(0.4476,1.821))
    mechanism2.omega_2 = omega_vec_out[3]
    #mechanism2.animate()
    #mechanism2.animate()

    print(numerical.numerical(th2,hy,3.07181611,3.99108686,3.05088475))
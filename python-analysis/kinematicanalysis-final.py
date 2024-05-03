from matplotlib import pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation
import numpy as np
import csv
import pandas as pd
from matplotlib.patches import Polygon

#import numerical

ground = np.hypot(1.821, 0.4476)
l1 = ground # ground 1
l2 = 3 # crank 1
l3 = 4 # coupler 1
l4 =  3 # ternary 1

l5 = 3 # ternary 2
l6 = 4 # coupler 2
l7 = 3 # end
l8 = ground

origin_y = 0.4476

def getBoundedAngle(angle):
    if angle < 0:
        return 360+angle
    elif angle > 360:
        return angle-360
    else:
        return angle

def fourbarpos(a, b, c, d, th2_array, delta = -1): # default for open configuration
    th3_arr = []
    th4_arr = []

    for i in th2_array:
    
        # K constants based on link lengths
        K1 = d / a
        K2 = d / c
        K3 = (a ** 2 - b ** 2 + c ** 2 + d ** 2) / (2 * a * c)
        K4 = d / b
        K5 = (c ** 2 - d ** 2 - a ** 2 - b ** 2) / (2 * a * b)
        
        th2 = np.deg2rad(i)  # converting to rads

        A = -K1 - K2 * np.cos(th2) + K3 + np.cos(th2)
        B = -2 * np.sin(th2)
        C = K1 - K2 * np.cos(th2) + K3 - np.cos(th2)
        D = np.cos(th2) - K1 + K4 * np.cos(th2) + K5
        E = -2 * np.sin(th2)
        F = K1 + (K4 - 1) * np.cos(th2) + K5
        
        # We don't want to handle imaginary roots
        disc_4 = B ** 2 - (4 * A * C)
        disc_3 = E ** 2 - (4 * D * F)
        if disc_4 < 0 or disc_3 < 0:
            print (B, A, C, E, D, F)
            print (disc_3, disc_4)
            print('rip2')
            raise SystemExit('Error: This function does not handle imaginary roots')
        
        # Solve for thetas 
        th4 = 2 * np.arctan2((-B + delta * np.sqrt(B ** 2 - 4 * A * C)), (2 * A))
        th3 = 2 * np.arctan2((-E + delta * np.sqrt(E ** 2 - 4 * D * F)), (2 * D))

        th3_arr.append(getBoundedAngle(np.rad2deg(th3)))
        th4_arr.append(getBoundedAngle(np.rad2deg(th4)))

    return np.array(th3_arr), np.array(th4_arr)
    
def fourbarvel(a, b, c, d, th2_arr, omega2_arr):

    th3_arr, th4_arr = fourbarpos(a, b, c, d, th2_arr)

    # initialize angular velocities
    omega3_arr = np.array([0.0 for _ in range(len(th2_arr))])
    omega4_arr = np.array([0.0 for _ in range(len(th2_arr))]) 

    # initialize linear velocities
    VA = np.zeros((len(th2_arr), 2))
    VBA = np.zeros((len(th2_arr), 2))
    VB = np.zeros((len(th2_arr), 2))
  
    for i in range(len(th2_arr)):
      omega2 = omega2_arr[i]
      th2 = np.deg2rad(th2_arr[i])
      th3 = np.deg2rad(th3_arr[i])
      th4 = np.deg2rad(th4_arr[i])

      omega3 = omega2 * (a/b) * (np.sin(th4 - th2)) / (np.sin(th3 - th4))
      omega4 = omega2 * (a/c) * (np.sin(th2 - th3)) / (np.sin(th4 - th3))

      # linear velocities of the joints
      VA[i] = np.array([a * omega2 * x for x in [-np.sin(th2),np.cos(th2)]])
      VBA[i] = np.array([b * omega3 * x for x in [-np.sin(th3),np.cos(th3)]])
      VB[i] = np.array([c * omega4 * x for x in [-np.sin(th4),np.cos(th4)]])

      omega3_arr[i] = omega3
      omega4_arr[i] = omega4

    return omega3_arr, omega4_arr, VA, VBA, VB

def fourbaraccel(a, b, c, d, th2_arr, omega2_arr, alpha2_arr):

    th3_arr, th4_arr = fourbarpos(a, b, c, d, th2_arr)
    omega3_arr, omega4_arr, VA, VBA, VB = fourbarvel(a, b, c, d, th2_arr, omega2_arr)

    alpha3_arr = np.array([0.0 for _ in range(len(th2_arr))])
    alpha4_arr = np.array([0.0 for _ in range(len(th2_arr))])

    AA = np.zeros((len(th2_arr), 2))
    ABA = np.zeros((len(th2_arr), 2))
    AB = np.zeros((len(th2_arr), 2))

    for i in range(len(th2_arr)):
        th2 = np.deg2rad(th2_arr[i])
        th3 = np.deg2rad(th3_arr[i])
        th4 = np.deg2rad(th4_arr[i])
        omega2 = omega2_arr[i]
        omega3 = omega3_arr[i]
        omega4 = omega4_arr[i]
        alpha2 = alpha2_arr[i]

        A = c * np.sin(th4)
        B = b * np.sin(th3)
        C = a * alpha2 * np.sin(th2) + a * omega2 ** 2 * np.cos(th2) + b * omega3 ** 2 * np.cos(th3) - c * omega4 ** 2 * np.cos(th4)

        D = c * np.cos(th4)
        E = b * np.cos(th3)
        F = a * alpha2 * np.cos(th2) + a * omega2 ** 2 * np.sin(th2) + b * omega3 ** 2 * np.sin(th3) - c * omega4 ** 2 * np.sin(th4)

        alpha3 = (C * D - A * F) / (A * E - B * D)
        alpha4 = (C * E - B * F) / (A * E - B * D)

        A_Ax = -a * alpha2 * np.sin(th2) - a * omega2 ** 2 * np.cos(th2)
        A_Ay = a * alpha2 * np.cos(th2) - a * omega2 ** 2 * np.sin(th2)
        A_BAx = -b * alpha3 * np.sin(th3) - b * omega3 ** 2 * np.cos(th3)
        A_BAy = b * alpha3 * np.cos(th3) - b * omega3 ** 2 * np.sin(th3)
        A_Bx = -c * alpha4 * np.sin(th4) - c * omega4 ** 2 * np.cos(th4)
        A_By = c * alpha4 * np.cos(th4) - c * omega4 ** 2 * np.sin(th4)

        alpha3_arr[i] = alpha3
        alpha4_arr[i] = alpha4
        AA[i] = [A_Ax, A_Ay]
        ABA[i] = [A_BAx, A_BAy]
        AB[i] = [A_Bx, A_By]

    return alpha3_arr, alpha4_arr, AA, ABA, AB

def point(a, th2_arr, th3_arr, omega3_arr, VA_arr, y = 2.55, x = 1.8):
    delta = np.arctan2(x, y)
    p = np.hypot(x, y)
    Vpa = np.zeros((len(th2_arr), 2))
    Vp = np.zeros((len(th2_arr), 2))
    p_arr = np.zeros((len(th2_arr), 2))
    for i in range(len(th3_arr)):
        th3 = np.deg2rad(th3_arr[i])
        th2 = np.deg2rad(th2_arr[i])
        omega3 = omega3_arr[i]
        VA = VA_arr[i]

        Px = a * np.cos(th2) + p * np.cos(th3 + delta)
        Py = a * np.sin(th2) + p * np.sin(th3 + delta)
        print(Px,Py)
        p_arr[i] = np.array([Px,Py])

        Vpa[i] = np.array([p * omega3 * x for x in [-np.sin(th3 + delta), np.cos(th3 + delta)]])
        Vp[i] = np.array([Vpa[i][0] + VA[0], Vpa[i][1] + VA[1]])
    return p_arr, Vp

th2 = np.linspace(0, 360, 360)
#th2 -= 13.8
th3, th4 = (fourbarpos(l2, l3, l4, l1, th2))
input_angle = th4 - 38.9 + 13.8
th5, th6 = (fourbarpos(l5, l6, l7, l8, input_angle))
omega2 = np.ones_like(th2) * -1
alpha2 = np.zeros_like(th2)

omega3, omega4, VA, VBA, VB = fourbarvel(l2, l3, l4, l1, th2, omega2)
omega5, omega6, VC, VDC, VD = fourbarvel(l5,l6,l7,l8, input_angle, omega4)

alpha3, alpha4, AA, ABA, AB = fourbaraccel(l2, l3, l4, l1, th2, omega2, alpha2)
alpha5, alpha6, AC, ADC, AD = fourbaraccel(l5, l6, l7, l8, input_angle, omega4, alpha4)

p, Vp = point(l5, input_angle, th5, omega5, VC, 2.55)
vel = [np.hypot(v[0], v[1]) for v in Vp]

pb, Vbp = point(l5, th4 - 38.9, th5, omega5, VC, y = 0)
velb = [np.hypot(v[0], v[1]) for v in Vp]

def plots():
    plt.plot(th2, omega6)
    plt.xlabel('Theta 2')
    plt.ylabel('Angular velocity of theta 6')
    plt.title('θ2 vs ω6')
    plt.savefig('omega6th2.png')
    plt.show()

    rin = l2
    rout = l7

    w = omega2/omega6
    r = rin/rout

    plt.plot(th2, w*r)
    plt.xlabel('Theta 2')
    plt.title('Mechanical Advantage')
    plt.savefig('mechadv.png')
    plt.show()

    plt.plot(th2, alpha6)
    plt.xlabel('Theta 2')
    plt.ylabel('Angular Acceleration of Link 6')
    plt.title('Theta 2 vs Angular Acceleration of Link 6')
    plt.savefig('th2accel6.png')
    plt.show()

    plt.plot(th2, omega6/omega2)
    plt.xlabel('Theta 2')
    plt.ylabel('Angular Velocity Ratio')
    plt.title('Theta 2 vs Velocity Ratio')
    plt.savefig('th2velratio.png')
    plt.show()

    for i in np.linspace(0, 3, 6):
        p, Vp = point(l5, th4 - 39.01, th5, omega5, VC, i)
        Px = p[:, 0]
        Py = p[:, 1]
        plt.plot(Px, Py)
    #Px = p[:, 0]
    #Py = p[:, 1]
    #plt.plot(Px, Py)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('End Effector Curve')
    plt.legend()
    plt.savefig('curve2.png')
    plt.show()

    plt.plot(th2, vel)
    plt.xlabel('Theta 2')
    plt.ylabel('Velocity of end effector [in/s]')
    plt.title('End Effector Velocity')
    plt.savefig('endvel2.png')
    plt.show()
plots()

def save():
    csv_file_path = "solution.csv"

    new_data = pd.DataFrame({
        'Theta2': th2,
        'Theta3': th3,
        'Theta4': th4,
        'Theta5': th5,
        'Theta6': th6,
        'Omega2': omega2,
        'Omega3': omega3,
        'Omega4': omega4,
        'Omega5': omega5,
        'Omega6': omega6,
        'Alpha2': alpha2,
        'Alpha3': alpha3,
        'Alpha4': alpha4,
        'Alpha5': alpha5,
        'Alpha6': alpha6,
    })
    new_data.to_csv(csv_file_path, mode='a', header=False, index=False)
'''
# https://www.youtube.com/watch?v=nT16-yQrnFk
fig, axis = plt.subplots()
animated_plt, = axis.plot([], [])

axis.set_xlim([-5, 5])
axis.set_ylim([-5, 5])
axis.set_title('Animation of Golf Mechanism')

plt.grid()

t = np.linspace(0, 360, 361)

def update(frame):
    animated_plt.set_data()
    return

animation = FuncAnimation(fig = fig, 
                          func = update,
                          frames = len(t)
                          interval=)'''

# Make sure numpy and matplotlib are correctly imported
# and that you have defined all necessary link lengths and angles arrays

def calculate_positions(l2, l3, l4, l5, l6, l7, l8, th2, th3, th4, th5, th6, origin):
    # Convert angles to radians
    ang = np.deg2rad(13.8)
    th2_rad = np.deg2rad(th2)
    th3_rad = np.deg2rad(th3)
    th4_rad = np.deg2rad(th4) 
    th4_2rad = np.deg2rad(th4 - 38.9 - 13.8) 
    th5_rad = np.deg2rad(th5)
    th6_rad = np.deg2rad(th6)
    

    O = [-origin[0],                                0]
    A = [-origin[0] + l2 * np.cos(th2_rad),         l2 * np.sin(th2_rad)]
    B = [l4 * np.cos(th4_rad),                      origin[1] + l4 * np.sin(th4_rad)] # grounded
    C = [0,                                         origin[1]]
    D = [l5 * np.cos(th4_2rad),                     origin[1] + l5 * np.sin(th4_2rad)] # grounded
    E = [origin[0] + l7 * np.cos(th6_rad),          l7 * np.sin(th6_rad)]
    F = [origin[0],                                 0]

    return [O, A, B, C, D, E, F]

fig, ax = plt.subplots()
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_aspect('equal')

origin = [1.821, 0.4476]  # Central point 
# Initialize th2, th3, th4, th5, th6 here as needed

lines = [
    ax.plot([], [], 'o-', lw=2, color='red')[0],    # Line between O and A
    ax.plot([], [], 'o-', lw=2, color='yellow')[0],   # Line between A and B
    ax.plot([], [], 'o-', lw=2, color='blue')[0],  # Line between B and C
    ax.plot([], [], 'o-', lw=2, color='blue')[0], # Line between B and D
    ax.plot([], [], 'o-', lw=2, color='blue')[0], # Line between C and D
    ax.plot([], [], 'o-', lw=2, color='orange')[0], # Line between D and E
    ax.plot([], [], 'o-', lw=2, color='cyan')[0],   # Line between E and F
    ax.plot([], [], 'o-', lw=2, color='black')[0],   # Line between O and C
    ax.plot([], [], 'o-', lw=2, color='black')[0],  # Line between C and F
]


def update(frame):
    positions = calculate_positions(l2, l3, l4, l5, l6, l7, l8, th2[frame], th3[frame], th4[frame], th5[frame], th6[frame], origin)
    O, A, B, C, D, E, F = positions
    # Update line data
    lines[0].set_data([O[0], A[0]], [O[1], A[1]])  # O-A
    lines[1].set_data([A[0], B[0]], [A[1], B[1]])  # A-B
    lines[2].set_data([B[0], C[0]], [B[1], C[1]])  # B-C
    lines[3].set_data([B[0], D[0]], [B[1], D[1]])  # B-D
    lines[4].set_data([C[0], D[0]], [C[1], D[1]])  # C-D
    lines[5].set_data([D[0], E[0]], [D[1], E[1]])  # D-E
    lines[6].set_data([E[0], F[0]], [E[1], F[1]])  # E-F
    lines[7].set_data([O[0], C[0]], [O[1], C[1]])  # O-C
    lines[8].set_data([C[0], F[0]], [C[1], F[1]])  # C-F

    # Update the polygon to shade the area formed by B, C, D
    polygon.set_xy([B, C, D])
    return lines + [polygon]

lines = [ax.plot([], [], 'o-', lw=2)[0] for _ in range(9)]
polygon = Polygon([[0, 0], [0, 0], [0, 0]], closed=True, color='blue', alpha=0.5)
ax.add_patch(polygon)

def init():
    for line in lines:
        line.set_data([], [])
    polygon.set_xy([[0, 0], [0, 0], [0, 0]])
    return lines + [polygon]

ani = FuncAnimation(fig, update, frames=np.arange(len(th2)), init_func=init, blit=True)
plt.show()

ani.save('anim.gif', dpi=200, writer='pillow')

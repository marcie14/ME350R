import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time

global ani

# given
ground = np.hypot(1.821, 0.4476)
l1 = ground # ground 1
l2 = 3.04115625 # crank 1
l3 = 4 # coupler 1
l4 =  3.06065793 # ternary 1

l5 = 3.06065793 # ternary 2
l6 = 4 # coupler 2
l7 = 3.04115625 # end
l8 = ground

# ground geometry (for plotting)
verts = [[0,0], [3.64, 0], [1.821, 0.4476], [0,0]]

origin_y = 0.4476

omega_2 = -1

def fourbarpos(a,b,c,d,th_1,th2,delta=-1):
    th_3_arr = []
    th_4_arr = []
    # print(th2)

    for i in (th2):
    
        # K constants based on link lengths
        K1 = d/a
        K2 = d/c 
        K3 = (a ** 2 - b ** 2 + c ** 2 + d ** 2) / (2 * a * c)
        K4 = d/b
        K5 = (c ** 2 - d ** 2 - a ** 2 - b ** 2)/(2 * a * b)
        
        th_1 = np.deg2rad(th_1)
        th_2 = np.deg2rad(i)  # converting to rads
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
            # print (B,A,C,E,D,F)
            # print (disc_3,disc_4)
            # print('rip2')
            raise SystemExit('Error: This function does not handle imaginary roots')
        
        # Solve for thetas 
        th_4 = 2*np.arctan2((-B + delta*np.sqrt(B**2-4*A*C)),(2*A))
        th_3 = 2*np.arctan2((-E + delta*np.sqrt(E**2-4*D*F)),(2*D))

        th_3_arr.append(np.rad2deg(th_3))
        th_4_arr.append(np.rad2deg(th_4))

    return np.zeros(361), th2, np.array(th_3_arr), np.array(th_4_arr)

def plotFourBar(links, thetas, index):
    assert len(links) == 4, "There should be 4 link lengths"
    assert len(thetas) == 4, "There should be 4 lists of angles"

    # Convert angles to radians
    thetas = [np.deg2rad(theta[index]) for theta in thetas]

    # Calculate positions
    x_positions = [0]
    y_positions = [0] 

    for i in range(2):
        x_positions.append(x_positions[i] + links[i]*np.cos(thetas[i]))
        y_positions.append(y_positions[i] + links[i]*np.sin(thetas[i]))

    # # For the last link, we need to close the loop
    x_positions.append(1.821)
    y_positions.append(0.4476)
    
    x_positions.append(0)
    y_positions.append(0)


    # print(x_positions)
    # print(y_positions)
    
    # Plot
    plt.figure()
    
    # plot ground triangle
    gndX, gndY = zip(*verts)
    plt.plot(gndX, gndY, 'r-', lw=2)
    plt.fill(gndX, gndY, 'r', alpha=0.3)
    
    for i in range(4):
        plt.plot([x_positions[i], x_positions[i+1]], [y_positions[i], y_positions[i+1]], 'o-', lw=2)
    
    plt.xlim(-4, 8)
    plt.ylim(-6, 6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.legend(['Ground', 'Ground', 'Link 2', 'Link 3', 'Link 4', 'Link 1'])
    plt.grid()
    plt.show()
    
def animateFourBar(links, thetas):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    # Initial positions
    x_positions = [0]
    y_positions = [0]

    for i in range(2):
        x_positions.append(x_positions[i] + links[i]*np.cos(thetas[i]))
        y_positions.append(y_positions[i] + links[i]*np.sin(thetas[i]))

    # # For the last link, we need to close the loop
    # x_positions.append(1.821)
    # y_positions.append(0.4476)

    # Set the x and y limits
    ax.set_xlim(-4, 8)
    ax.set_ylim(-6, 6)

    # Plot
    line, = ax.plot(x_positions[0], y_positions[0], 'o-', lw=2)

    def update(frame):
        # Update thetas
        thetas_frame = [np.deg2rad(theta[frame]) for theta in thetas]

        # Calculate new positions
        x_positions = [0] 
        y_positions = [0] 
        
        # plot ground triangle
        gndX, gndY = zip(*verts)
        plt.plot(gndX, gndY, 'r-', lw=2)
        plt.fill(gndX, gndY, 'r', alpha=0.3)
        
        plt.grid()
        for i in range(2):
            x_positions.append(x_positions[i] + links[i]*np.cos(thetas_frame[i]))
            y_positions.append(y_positions[i] + links[i]*np.sin(thetas_frame[i]))


        # For the last link, we need to close the loop
        x_positions.append(1.821)
        y_positions.append(0.4476)
    
        x_positions.append(0)
        y_positions.append(0)
        
        # print(x_positions[2], y_positions[2])
        # time.sleep(0.2)
        # Update line data
        line.set_data(x_positions, y_positions)
        return line,
    
    ani = FuncAnimation(fig, update, frames=range(len(thetas[0])), blit=True, interval = 10)

    plt.show()

def plotSixBar(links, thetas, index):
    assert len(links) == 8, "There should be 8 link lengths"
    assert len(thetas) == 8, "There should be 8 lists of angles"

    # Convert angles to radians
    thetas = [np.deg2rad(theta[index]) for theta in thetas]

    # Calculate positions
    x_positions = [0]
    y_positions = [0] 

    for i in range(2):
        x_positions.append(x_positions[i-1] + links[i]*np.cos(thetas[i]))
        y_positions.append(y_positions[i-1] + links[i]*np.sin(thetas[i]))

    # # For the last link, we need to close the loop
    x_positions.append(1.821)
    y_positions.append(0.4476)
    
    # x_positions.append(0)
    # y_positions.append(0)
    for i in range(4, 6):
        x_positions.append(x_positions[i-1] + links[i]*np.cos(thetas[i]))
        y_positions.append(y_positions[i-1] + links[i]*np.sin(thetas[i]))
    
    # print(len(x_positions))

    x_positions.append(3.64293832)
    y_positions.append(0)

    print(x_positions)
    print(y_positions)
    
    # Plot
    plt.figure()
    
    # plot ground triangle
    gndX, gndY = zip(*verts)
    plt.plot(gndX, gndY, 'r-', lw=2)
    plt.fill(gndX, gndY, 'r', alpha=0.3)
    print(len(x_positions))
    for i in range(6):
        plt.plot([x_positions[i], x_positions[i+1]], [y_positions[i], y_positions[i+1]], 'o-', lw=2)
    
    plt.xlim(-4, 8)
    plt.ylim(-6, 6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.legend(['Ground', 'Ground', 'Link 2', 'Link 3', 'Link 4', 'Link 4', 'Link 5', 'Link 6'])
    plt.grid()
    plt.show()
    
def animateSixBar(links, thetas):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    # Initial positions
    x_positions = [0]
    y_positions = [0]

    for i in range(2):
        for j in range(2):
            x_positions.append(x_positions[i] + links[i]*np.cos(thetas[i][j]))
            y_positions.append(y_positions[i] + links[i]*np.sin(thetas[i][j]))

    # # For the last link, we need to close the loop
    # x_positions.append(1.821)
    # y_positions.append(0.4476)

    # Set the x and y limits
    ax.set_xlim(-4, 8)
    ax.set_ylim(-6, 6)
    # print(1)
    # Plot
    line, = ax.plot(x_positions[0], y_positions[0], 'o-', lw=2)
    # print(type(thetas))
    def update(frame):
        # Update thetas
        # print(frame)

        thetas_frame = [np.deg2rad(theta[frame]) for theta in thetas]
        thetas_frame = [list(x) for x in zip(*thetas_frame)]
        print(thetas_frame[0][0])

        # Calculate new positions
        x_positions = [0] 
        y_positions = [0] 
        
        # plot ground triangle
        gndX, gndY = zip(*verts)
        plt.plot(gndX, gndY, 'r-', lw=2)
        plt.fill(gndX, gndY, 'r', alpha=0.3)
        # print(2)
        plt.grid()
        for i in range(2):
            x_positions.append(x_positions[i] + links[i]*np.cos(thetas_frame[frame][i]))
            y_positions.append(y_positions[i] + links[i]*np.sin(thetas_frame[frame][i]))

        # print(3)
        # For the last link, we need to close the loop
        x_positions.append(1.821)
        y_positions.append(0.4476)
    
        # x_positions.append(0)
        # y_positions.append(0)
        # print('x pos ' + str(len(x_positions)))
        # print('test', links[4])
        for i in range(4, 6):
            
            x_positions.append(x_positions[i-1] + links[i]*np.cos(thetas_frame[frame][i]))
            y_positions.append(y_positions[i-1] + links[i]*np.sin(thetas_frame[frame][i]))
    
    
        # print('x pos ' + str(len(x_positions)))
        # print(x_positions[3], y_positions[3], np.hypot(x_positions[3], y_positions[3]))
        # print(len(x_positions))

        x_positions.append(3.64293832)
        y_positions.append(0)
        # print(4)
        
        # print(x_positions[2], y_positions[2])
        # time.sleep(0.2)
        # Update line data
        line.set_data(x_positions, y_positions)
        return line,
    
    ani = FuncAnimation(fig, update, frames=range(len(thetas[0])), blit=True, interval = 10)

    plt.show()



def gptanimateSixBar(links, thetas, velocities):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    # Initial positions
    x_positions = [0]
    y_positions = [0]

    for i in range(2):
        for j in range(2):
            x_positions.append(x_positions[i] + links[i]*np.cos(thetas[i][j]))
            y_positions.append(y_positions[i] + links[i]*np.sin(thetas[i][j]))

    ax.set_xlim(-4, 8)
    ax.set_ylim(-4, 6)
    Q = ax.quiver([], [], [], [], angles='xy', scale_units='xy', scale=1)

    line, = ax.plot(x_positions[0], y_positions[0], 'o-', lw=2)

    def update(frame):
        # thetas_frame = [np.deg2rad(theta[frame]) for theta in thetas]
        thetas_frame = [np.deg2rad(theta[frame]) for theta in thetas]
        velocity_frame = [vel[frame] for vel in velocities]
        print(links)
        # thetas_frame = [list(x) for x in zip(*thetas_frame)]

        x_positions = [0] 
        y_positions = [0] 

        gndX, gndY = zip(*verts)
        plt.plot(gndX, gndY, 'r-', lw=2)
        plt.fill(gndX, gndY, 'r', alpha=0.3)
        plt.grid()
        
        print(f"Frame: {frame}, Length of thetas[0]: {len(thetas[0])}")
        print(f"Length of thetas_frame: {len(thetas_frame)}")
        print(f"Length of links: {len(links)}")
        
        for i in range(2):
            x_positions.append(x_positions[i] + links[i]*np.cos(thetas_frame[i]))
            y_positions.append(y_positions[i] + links[i]*np.sin(thetas_frame[i]))

        # hardcode the ground link
        x_positions.append(1.821)
        y_positions.append(0.4476)

        for i in range(4, 6):
            x_positions.append(x_positions[i-1] + links[i]*np.cos(thetas_frame[i]))
            y_positions.append(y_positions[i-1] + links[i]*np.sin(thetas_frame[i]))

        # hardcode the ground link
        x_positions.append(3.64293832)
        y_positions.append(0)

        U = [vel[0] for vel in velocity_frame]
        V = [vel[1] for vel in velocity_frame]
        print(len(U), len(V))
        print(len(x_positions), len(y_positions))
        Q.set_UVC(U, V)
        Q.set_offsets(np.array([x_positions, y_positions]).T)
        # for i in range(len(x_positions) - 1):
        #     # Calculate the velocity components
        #     u = velocity_frame[i]
        #     v = velocity_frame[i]

        #     # Plot the arrow
        #     plt.quiver(x_positions[i], y_positions[i], u, v, angles='xy', scale_units='xy', scale=1)
        
        line.set_data(x_positions, y_positions)
        return line, Q,

    ani = FuncAnimation(fig, update, frames=range(len(thetas[0])), blit=True, interval = 10)

    plt.show()
    return ani  # Return the animation object



########## MAIN ##########
    
th1, th2, th3, th4 = (fourbarpos(l2,l3,l4,l1,0,np.linspace(0,361,361)))
th1 = [13.8094803531 for x in th1] # angle of ground
links1 = [l2, l3, l4, l1]
thetas1 = [th2, th3, th4, th1]
# plotFourBar(links1, thetas1, 113)
# animateFourBar(links1, thetas1)


th5, th6, th7, th8 = (fourbarpos(l5,l6,l7,l8,0,th4 - 39.01))
# print(th4 - th6) # th6 is constant angle between ternary links
links2 = [l5, l6, l7, l8] # note that position of l5 is not the same as th5 in list below
thetas2 = [th6, th7, th8, th5]
# plotFourBar(links2, thetas2, 113)
# animateFourBar(links2, thetas2)

links3 = [l2, l3, l4, l1, l5, l6, l7, l8] # note that position of l5 is not the same as th5 in list below
# print(links3)
thetas3 = [th2, th3, th4, th1, th6, th7, th8, th5]
# plotSixBar(links3, thetas3, 50)
# animateSixBar(links3, thetas3)

vel1 = [[1, 1] for x in range(len(th1))]
vel2 = [[1, 1] for x in range(len(th2))]
vel3 = [[1, 1] for x in range(len(th3))]
vel4 = [[1, 1] for x in range(len(th4))]
vel5 = [[1, 1] for x in range(len(th5))]
vel6 = [[1, 1] for x in range(len(th6))]
vel7 = [[1, 1] for x in range(len(th7))]

vels = [vel2, vel3, vel4, vel5, vel6, vel7, vel1]
gptanimateSixBar(links3, thetas3, vels)


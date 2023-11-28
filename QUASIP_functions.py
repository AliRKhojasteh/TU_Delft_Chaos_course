
#functions

import matplotlib.pyplot as plt
import numpy as np


# Define the main function
def quasip_test(Map, alpha, beta, gamma, A, B, Omega, Picture, k, b):
    global h
    animsteps = 32

    Map = 0 

    # Initialization of variables
    h = 2 * np.pi / (animsteps * Omega)  # Step size calculation
    # print('Step size is ', h)
    
    # Initial state vector
    # xc = np.array([0.3 * 2 * np.pi, 0.3, 0])
    xc = np.array([0.3 * 2 * np.pi, 0.3, 0])

    # yy = xc[0] / (2 * np.pi)  # Normalizing the first element of xc
    yy = [xc[0] / (2 * np.pi), 0]  # Normalizing the first element of xc and initializing the second element

    tn = yy[0]
    rn = 0.3
    t = 0

    # Warning if Map is on in combination with Picture 3
    if Map == 1 and Picture == 3:
        print('The option Map is not compatible with Picture option 3, Picture option is set to 1')
        Picture = 1
    
    n = 1
    nit = 20
    positions_1 = []
    positions_2 = []
    positions_3 = []
    positions_4 = []

    while n < nit:
        if Map == 0:
            b1 = np.sin(xc[0])
            b2 = -np.cos(xc[0])

            # print('b1 = ', b1, 'b2 = ', b2)
            to = tn
            # xx=0
            xx = [0, 0]  # Initializing xx as a list with two elements

            # Update state using Runge-Kutta method
            for i in range(animsteps):
                t += h
                xc = Runge_test(xc,alpha,beta,gamma,A,B,Omega,h)
                b1 = np.sin(xc[0])
                b2 = -np.cos(xc[0])            
                positions_1.append((b1, b2))  # Append a tuple of (b1, b2) to the list

                
                if Picture == 3:
                    # 2D torus plot (theta vs time)
                    tn = xc[0] / (2 * np.pi)
                    yy[1] = tn % 1
                    xx[1] = xx[0] + h * Omega / (2 * np.pi)

                    if yy[1] < 0:
                        yy[1] += 1.0

                    # Periodic boundary conditions
                    if yy[1] > (yy[0] + 0.5):
                        xtemp = xx[0] - yy[0] * (xx[1] - xx[0]) / (yy[1] - 1 - yy[0])
                        xline = [xx[0], xtemp]
                        yline = [yy[0], 0]
                        # axs[1].plot(xline, yline, '-k')
                        positions_3.append((xline, yline))
                        xx[0], yy[0] = xtemp, 1
                    elif yy[1] < (yy[0] - 0.5):
                        xtemp = xx[0] + (1 - yy[0]) * (xx[1] - xx[0]) / (yy[1] + 1 - yy[0])
                        xline = [xx[0], xtemp]
                        yline = [yy[0], 1]
                        # axs[1].plot(xline, yline, '-k')
                        positions_3.append((xline, yline))
                        xx[0], yy[0] = xtemp, 0

                    # Update the lines for the current frame
                    xline = [xx[0], xx[1]]
                    yline = [yy[0], yy[1]]
                    # axs[1].plot(xline, yline, '-k')
                    positions_3.append((xline, yline))
                    xx[0], yy[0] = xx[1], yy[1]


            tn = xc[0] / (2 * np.pi)
            rn = xc[1]
        else:
            # Map is switched on
            ro = rn
            to = tn
            rn = b * ro - k / (2 * np.pi) * np.sin(2 * np.pi * (to % 1))
            tn = to + Omega + rn
            
        
        # Check if Picture is set to 1
        if Picture == 1:
            yy[0] = tn % 1
            xx[0] = to % 1
            
            positions_2.append((xx[0], yy[0]))
            # # Create subplot only if Map is 0
            # if Map == 0:
            #     ax = plt.subplot(1, 2, 2)
            # else:
            #     ax = plt.gca()  # Get current axis

            # Plot the 1D map
            # ax.plot(xx, yy, '.k')
            # plt.pause(0.01)  # Pause to update the plot
        elif Picture == 2:
            yy[0] = rn
            xx[0] = tn % 1
            positions_4.append((xx[0], yy[0]))

        n += 1        

    return positions_1, positions_2, positions_4, positions_3






def equations(x,alpha,beta,gamma,A,B,Omega):
    # Define the differential equations
    f = np.zeros(3)
    f[0] = x[1]
    f[1] = 1 / alpha * (-beta * x[1] - gamma * np.sin(x[0]) + A + B * np.cos(Omega * x[2]))
    f[2] = 1
    return f

def Runge_test(xc,alpha,beta,gamma,A,B,Omega,h):
    # Runge-Kutta 4th order method
    n = len(xc)
    x = np.zeros_like(xc)
    c1 = np.zeros_like(xc)
    c2 = np.zeros_like(xc)
    c3 = np.zeros_like(xc)
    c4 = np.zeros_like(xc)

    for i in range(n):
        x[i] = xc[i]
    f = equations(x,alpha,beta,gamma,A,B,Omega)
    for i in range(n):
        c1[i] = h * f[i]

    for i in range(n):
        x[i] = xc[i] + c1[i] / 2
    f = equations(x,alpha,beta,gamma,A,B,Omega)
    for i in range(n):
        c2[i] = h * f[i]

    for i in range(n):
        x[i] = xc[i] + c2[i] / 2
    f = equations(x,alpha,beta,gamma,A,B,Omega)
    for i in range(n):
        c3[i] = h * f[i]

    for i in range(n):
        x[i] = xc[i] + c3[i]
    f = equations(x,alpha,beta,gamma,A,B,Omega)
    for i in range(n):
        c4[i] = h * f[i]

    for i in range(n):
        xc[i] = xc[i] + (c1[i] + 2 * c2[i] + 2 * c3[i] + c4[i]) / 6

    return xc

def stopf(event):
    # Event handler for the stop button
    # In a GUI application, this would terminate the event loop or close the window
    pass  # Placeholder, functionality depends on GUI framework used



def init_figure():
    # Clear current figure
    plt.clf()

    # Creating a 'Stop' button - this is a placeholder and might need further integration with GUI event loop
    # For full functionality, you might need to use a GUI framework like Tkinter
    stop_button = plt.axes([0.81, 0.05, 0.1, 0.075])
    button = plt.Button(stop_button, 'Stop')
    button.on_clicked(stopf)  # Linking to a stop function

    if Map == 0:
        # Create first subplot
        ax1 = plt.subplot(1, 2, 1)
        ax1.axis([-1.5, 1.5, -1.5, 1.5])
        ax1.set_aspect('equal', 'box')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_title('Pendanim')

        # Create second subplot
        ax2 = plt.subplot(1, 2, 2)
    else:
        ax2 = plt.subplot(1, 1, 1)

    ax2.axis([0, 1, 0, 1])
    ax2.set_aspect('equal', 'box')

    # Set labels based on Picture value
    if Picture == 1:
        ax2.set_xlabel(r'$\theta_n$', fontsize=12)
        ax2.set_ylabel(r'$\theta_{n+1}$', fontsize=12)
    elif Picture == 2:
        ax2.set_xlabel(r'$\theta_n$', fontsize=12)
        ax2.set_ylabel(r'$\dot{\theta_n}$', fontsize=12)
        if Map == 0:
            ax2.axis([0, 1, -3, 3])
        else:
            ax2.axis([0, 1, -0.5, 0.5])
    elif Picture == 3:
        ax2.set_xlabel('t mod (2π / Ω)')
        ax2.set_ylabel('θ mod 2π')

    plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)  # Adjust subplot to make room for button

    fig = plt.gcf()
    fig.set_size_inches(20, 10)  # Set the figure size to 20x10

    
    while n < nit:
        if Map == 0:
            b1 = np.sin(xc[0])
            b2 = -np.cos(xc[0])

            # Update state using Runge-Kutta method
            for i in range(animsteps):
                t += h
                xc = Runge(xc)

                # Clear and update the plot for pendulum
                plt.subplot(1, 2, 1).cla()
                plt.plot(0, 0, '+', markersize=10)
                plt.plot([0, b1], [0, b2], 'b-')
                plt.plot(b1, b2, 'r.', markersize=25)
                plt.pause(0.01)

                # Additional plotting based on Picture
                if Picture == 3:
                    # Implement the logic for 2D torus plot
                    pass  # Placeholder for the 2D torus plot logic

            tn = xc[0] / (2 * np.pi)
            rn = xc[1]

        else:
            # Logic when Map is switched on
            ro = rn
            to = tn
            rn = b * ro - k / (2 * np.pi) * np.sin(2 * np.pi * (to % 1))
            tn = to + Omega + rn

        # Plotting based on Picture
        if Picture == 1:
            yy = tn % 1
            xx = to % 1
            plt.subplot(1, 2, 2).plot(xx, yy, '.k')
            plt.pause(0.01)
        elif Picture == 2:
            yy = rn
            xx = tn % 1
            plt.subplot(1, 2, 2).plot(xx, yy, '.k')
            plt.pause(0.01)

        n += 1

    return fig, ax1 if Map == 0 else None, ax2


#functions

import matplotlib.pyplot as plt
import numpy as np


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

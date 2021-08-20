import numpy as np
from tools import Potential, Update
import matplotlib.pyplot as plt
from importlib import reload
from time import sleep
from tqdm import tqdm
import matplotlib.animation as animation
import pickle
import sys

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

d = load_obj('all_data_{}part_{}years_{}pc'.format(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])))
d_taille = d['taille']
d_pc = d['pc']
rho = d['density']
Ndim = int(sys.argv[4])

if Ndim == 2 :
    d = d['position'][:, :, :2]

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.set(xlim=(-2*d_taille,2*d_taille), ylim=(-2*d_taille,2*d_taille))
    point=ax1.scatter([], [], s = 0.1, color = 'black')

    label2 = ax1.text(1.9*d_taille, 1.6*d_taille, "",
                ha='right', va='center',
                fontsize=10)

    ax2 = fig.add_subplot(122)

    im = ax2.imshow(rho[0])
    ttl1 = ax2.text(.2, 1.05, r'Trajectories as function of time $t$ (in kyr)', transform = ax1.transAxes, va='center')
    ttl2 = ax2.text(.4, 1.05, r'Number of stars $N$', transform = ax2.transAxes, va='center')

    def animate(i) :
        point.set_offsets(d[i]/d_pc)
        label2.set_text(i*int(sys.argv[2])/1e3)
        im.set_array(rho[i])
        return point, label2, im

    inter = 3
    ani = animation.FuncAnimation(fig, animate, frames=range(d.shape[0]), blit=True, interval=inter, repeat=True)
    plt.show()

save = sys.argv[5]

if Ndim == 3 :
    inter = 3
    # IMPORTS
    import numpy as np
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d.axes3d as p3
    import matplotlib.animation as animation

    def animate_scatters(iteration, data, scatters):
        """
        Update the data held by the scatter plot and therefore animates it.
        Args:
            iteration (int): Current iteration of the animation
            data (list): List of the data positions at each iteration.
            scatters (list): List of all the scatters (One per element)

        Returns:
            list: List of scatters (One per element) with new coordinates
        """
        for i in range(data[0].shape[0]):
            scatters[i]._offsets3d = (data[iteration][i,0:1], data[iteration][i,1:2], data[iteration][i,2:])
        return scatters

    def main(data, save=False):
        """
        Creates the 3D figure and animates it with the input data.
        Args:
        data (list): List of the data positions at each iteration.
        save (bool): Whether to save the recording of the animation. (Default to False).
        """

        # Attaching 3D axis to the figure
        fig = plt.figure()
        ax = p3.Axes3D(fig)

        # Initialize scatters
        scatters = [ ax.scatter(data[0][i,0:1], data[0][i,1:2], data[0][i,2:], s = 0.5, color = 'black') for i in range(data[0].shape[0]) ]

        # Number of iterations
        iterations = len(data)

        # Setting the axes properties
        #ax.set_xlim3d([-50, 50])
        ax.set_xlabel('X')

        #ax.set_ylim3d([-50, 50])
        ax.set_ylabel('Y')

        #ax.set_zlim3d([-50, 50])
        ax.set_zlabel('Z')

        ax.set_title('3D Animated Scatter Example')

        # Provide starting angle for the view.
        ax.view_init(25, 10)

        ani = animation.FuncAnimation(fig, animate_scatters, iterations, fargs=(data, scatters),
                                       interval=50, blit=False, repeat=True)

        if save:
            writergif = animation.PillowWriter(fps=35)
            ani.save('{}stars.gif'.format(int(sys.argv[1])),writer=writergif)
        else :
            plt.show()
    main(d['position'], save=save)

import numpy as np
import tools
import matplotlib.pyplot as plt
from importlib import reload
from time import sleep
from tqdm import tqdm
import pickle
import sys

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

during = int(sys.argv[1]) # Define the simulation time
existing = int(sys.argv[2]) # Allow to load the right number of years simulations
N = int(sys.argv[3])
Nyears = int(sys.argv[4])
seuil = int(sys.argv[5])


d = load_obj('all_data_{}part_{}years_{}pc'.format(N, existing*Nyears, seuil))

#------------------------------------------
# Loading parameters
N = int(d['position'].shape[1])
Nyears = d['Nyears']
dt = d['dt']
G = 6.67e-11/dt**2

pc = d['pc'] #m
Msun = 1.989e30 #kg
taille = d['taille']
L = d['L']
m = d['m']
nb = d['nb']
rho = np.zeros(((len(d['time']), nb, nb)))
Potential = tools.Potential(G)
Update = tools.Update(G)
#------------------------------------------




def density(X, t, nb) :

    border = np.linspace(-1.5*taille*pc, 1.5*taille*pc, nb)
    dens = np.zeros((len(border), len(border)))

    for i, j in enumerate(border):
        for k, l in enumerate(border):
            if i != nb-1 :
                if k != nb-1 :
                    ind = (X[t, :, 0] > border[i]) & (X[t, :, 0] < border[i+1]) & (X[t, :, 1] > border[k]) & (X[t, :, 1] < border[k+1])
                    dens[i, k] = np.sum(ind)
    return dens

def pos_init(N, r_lim):
    posi = np.zeros((N, 3))
    k=0
    while k < N :
        X = np.random.uniform(-r_lim, r_lim, (1, 3))

        #print(X)
        r = np.sqrt(X[0, 0]**2 + X[0, 1]**2 + X[0, 2]**2)
        if r < r_lim :
            #print(k)
            posi[k] = X
            k+=1
    return posi

old_pos = d['position']
old_vit = d['velocity']
old_rho = d['density']

new_pos = np.zeros(((during*d['position'].shape[0], N, 3)))
new_vit = np.zeros(((during*d['time'].shape[0], N, 3)))
new_rho = np.zeros(((during*len(d['time']), nb, nb)))

new_pos[:d['position'].shape[0]] = old_pos
new_vit[:d['velocity'].shape[0]] = old_vit
new_rho[:d['density'].shape[0]] = old_rho
#new_tab_pos = np.delete(new_tab_pos, len(d['time']), axis = 0)
#new_tab_vit = np.delete(new_tab_vit, len(d['time']), axis = 0)

time_array = np.arange(0, during*len(d['time'])*dt, dt)
#new_time = np.concatenate((d['time'], time_array), axis = 0)
#new_time = np.delete(new_time, len(d['time']))

print("\n\n======================== \n Computing trajectories \n======================== \n \n")

for t in tqdm(range(len(d['time']), len(time_array))) :
    #print(t)
    if t != len(time_array)-1:
        #print('computing new positions!')
        for i in range(N):
            r, ind, ax, ay, az = Update.update_acc(new_pos, m, t-1, i, seuil)

            if t == len(d['time']):
                new_pos[t, i, 0], new_vit[t, i, 0] = Update.first_step(new_pos[t-1, i, 0], new_vit[t-1, i, 0], ax, dt)
                new_pos[t, i, 1], new_vit[t, i, 1] = Update.first_step(new_pos[t-1, i, 1], new_vit[t-1, i, 1], ay, dt)
                new_pos[t, i, 2], new_vit[t, i, 2] = Update.first_step(new_pos[t-1, i, 2], new_vit[t-1, i, 2], az, dt)

            else:
                new_pos[t, i, 0] = Update.verlet(new_pos[t-1, i, 0], new_pos[t-2, i, 0], ax, dt)
                new_pos[t, i, 1] = Update.verlet(new_pos[t-1, i, 1], new_pos[t-2, i, 1], ay, dt)
                new_pos[t, i, 2] = Update.verlet(new_pos[t-1, i, 2], new_pos[t-2, i, 2], az, dt)

                new_vit[t, i, 0] = (new_pos[t-1, i, 0] - new_pos[t-2, i, 0])/dt
                new_vit[t, i, 1] = (new_pos[t-1, i, 1] - new_pos[t-2, i, 1])/dt
                new_vit[t, i, 2] = (new_pos[t-1, i, 2] - new_pos[t-2, i, 2])/dt

    # Last time step
    elif t == len(time_array)-1 :
        #print("last")
        for i in range(N):
            new_pos[t, i, 0] = new_pos[t-1, i, 0]
            new_pos[t, i, 1] = new_pos[t-1, i, 1]
            new_pos[t, i, 2] = new_pos[t-1, i, 2]

            new_vit[t, i, 0] = new_vit[t-1, i, 0]
            new_vit[t, i, 1] = new_vit[t-1, i, 1]
            new_vit[t, i, 2] = new_vit[t-1, i, 2]
    else:
        pass
        #print('No computing!')

    #old_pos = np.concatenate((old_pos, new_pos), axis = 0)
    #old_vit = np.concatenate((old_vit, new_vit), axis = 0)

print("\n\n======================== \n Computing density \n======================== \n \n")

for t in tqdm(range(len(d['time']), len(time_array))):
    if t >= len(d['time']):
        new_rho[t] = density(new_pos, t, nb)

new_dict = {
    'position' : new_pos,
    'velocity' : new_vit,
    'density' : new_rho,
    'taille' : taille,
    'L' : L,
    'pc' : pc,
    'Nyears' : Nyears,
    'dt' : dt,
    'seuil' : seuil,
    'nb' : nb,
    'time' : time_array,
    'm' : m}


def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

save_obj(new_dict, 'all_data_{}part_{}years_{}pc'.format(N, during*existing*Nyears, seuil))

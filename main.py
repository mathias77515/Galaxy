import numpy as np
import tools
import matplotlib.pyplot as plt
from importlib import reload
from time import sleep
from tqdm import tqdm
import pickle
import sys

Ndim = 3

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
    posi = np.zeros((N, Ndim))
    k=0
    while k < N :
        X = np.random.uniform(-r_lim, r_lim, (1, Ndim))
        #print(X)
        if Ndim == 3 :
            r = np.sqrt(X[0, 0]**2 + X[0, 1]**2 + X[0, 2]**2)
        elif Ndim == 2 :
            r = np.sqrt(X[0, 0]**2 + X[0, 1]**2)
        else:
            raise TypeError("Choose a 2D or 3D computation !")
        if r < r_lim :
            #print(k)
            posi[k] = X
            k+=1
    return posi

#------------------------------------------
N = int(sys.argv[1])
Nyears = int(sys.argv[4])
dt = 60*60*24*365*Nyears
G = 6.67e-11/dt**2

pc = 3.086e16 #m
Msun = 1.989e30 #kg
taille = 1e1
L = taille*pc/10

seuil_pc = float(sys.argv[3]) #pc
seuil = seuil_pc*pc

time = np.arange(0, int(sys.argv[2])*dt, dt)
nb = 60
rho = np.zeros(((len(time), nb, nb)))
Potential = tools.Potential(G)
Update = tools.Update(G)
#------------------------------------------

pos = np.zeros(((len(time), N, Ndim)))
p = pos_init(N, taille*pc)
pos[0] = p

m = np.random.uniform(0.7*Msun, 50*Msun,(N))
a = np.zeros(((len(time), N, Ndim)))
vit = np.zeros(((len(time), N, Ndim)))
vit[0] = np.random.normal(0, 1e3, (N, Ndim))


#vit[0] = np.zeros((N, 2))
r = np.zeros(((len(time), N, N)))
tab_fx = np.zeros(((len(time), N, N)))
tab_fy = np.zeros(((len(time), N, N)))
ind_all = np.arange(0, N)
i_f = 1

print("\n\n======================== \n Computing trajectories \n======================== \n \n")

for t in tqdm(range(len(time)-1)) :
    for i in range(N):
        r, ind, ax, ay, az = Update.update_acc(pos, m, t, i, seuil)

        if t == 0 :
            pos[t+1, i, 0], vit[t+1, i, 0] = Update.first_step(pos[t, i, 0], vit[t, i, 0], ax, dt)
            pos[t+1, i, 1], vit[t+1, i, 1] = Update.first_step(pos[t, i, 1], vit[t, i, 1], ay, dt)
            pos[t+1, i, 2], vit[t+1, i, 2] = Update.first_step(pos[t, i, 2], vit[t, i, 2], az, dt)
            #r_for_Ep = np.delete(r, i)
            #Ep_f[t+1] = np.sum(Potential.V(m_for_Ep, m[i], r_for_Ep))

        else:

            pos[t+1, i, 0] = Update.verlet(pos[t, i, 0], pos[t-1, i, 0], ax, dt)
            pos[t+1, i, 1] = Update.verlet(pos[t, i, 1], pos[t-1, i, 1], ay, dt)
            pos[t+1, i, 2] = Update.verlet(pos[t, i, 2], pos[t-1, i, 2], az, dt)

            vit[t+1, i, 0] = (pos[t+1, i, 0] - pos[t, i, 0])/dt
            vit[t+1, i, 1] = (pos[t+1, i, 1] - pos[t, i, 1])/dt
            vit[t+1, i, 2] = (pos[t+1, i, 2] - pos[t, i, 2])/dt



print("\n\n======================== \n Computing density \n======================== \n \n")

for t in tqdm(range(len(time)-1)):
    rho[t] = density(pos[:, :, :2], t, nb)

dict = {
    'position' : pos,
    'velocity' : vit,
    'density' : rho,
    'position' : pos,
    'taille' : taille,
    'L' : L,
    'pc' : pc,
    'Nyears' : Nyears,
    'dt' : dt,
    'seuil' : seuil,
    'nb' : nb,
    'time' : time,
    'm' : m}

def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

save_obj(dict, 'all_data_{}part_{}years_{}pc'.format(N, Nyears, int(seuil_pc)))

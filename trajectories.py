import multiprocessing as mp
import pandas as pnd
import matplotlib.pyplot as pt
import pytraj as pt
import os
import pandas as pnd
import json
import sys


def process_file(filename: str):
    # TODO: remove additional "and False" condition to activate this branch
    if filename[0][-16:-12] == "3EAK":
        pretrajectory_fn = folder_path + '/' + filename[0][-16:] + '/' + "mdcrd"
        trajectory_fn = folder_path + '/' + filename[0][-16:] + '/' + "amd.nc"
        prmtop_fn = folder_path + '/' + filename[0][-16:] + '/' + "ab.prmtop"
        print(filename[0][-12:])
        data = pnd.DataFrame({
            'mutant': filename[0][-11:],
            'total_radgyr': [],
            'total_rmsd': [],
            'molsurf': [],
            'CDR3_rmsd': [],
            'CDR3_radgyr': [],
            'CDR3_core_dist': [],
            'CDR3_40_dist': [],
            'CDR3_45_dist': [],
            '45_core_dist': []
        })
        traj = pt.superpose(pt.strip(pt.load([pretrajectory_fn, trajectory_fn], top=prmtop_fn, stride=10), ":WAT"))
        data['total_radgyr'] = pt.radgyr(traj, mask='@CA')
        data['CDR3_radgyr'] = pt.radgyr(traj, mask='103-118')
        data['total_rmsd'] = pt.rmsd(traj, mask='!@H')
        data['CDR3_rmsd'] = pt.rmsd(traj, mask="103-118")
        data['molsurf'] = pt.molsurf(traj, mask="!@H")
        data['CDR3_core_dist'] = pt.distance(traj, ":1-102,119-127 :103-118")
        data['CDR3_40_dist'] = pt.distance(traj, ":40 :103-118")
        data['CDR3_45_dist'] = pt.distance(traj, ":45 :103-118")
        data['45_core_dist'] = pt.distance(traj, ":1-102,119-127 :45")

        lock.acquire()
        try:
            total_data.append(data)
        finally:
            lock.release()


folder_path = '/home/misha/Dima/cameloid_AB/3EAK-derivatives/results'
container = os.walk(folder_path)
folder = []
for i in container:
    folder.append(i)

if __name__ == '__main__':
    total_data = []
    print('folder length = ' + str(len(folder)))
    mp.set_start_method('spawn')
    q = mp.Queue()
    p = mp.Process(target=process_file, args=(q,))
    p.start()
    print(q.get())
    p.join()
    print('before to pickl')
    pt.to_pickle(total_data, '/home/misha/Dima/cameloid_AB/3EAK-derivatives/results/total_data.pk')
    print(len(total_data))
    print(len(folder))

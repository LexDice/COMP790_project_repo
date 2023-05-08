# import trimesh
import os
import glob
import numpy as np
import pyvista as pv


dataset_path = '/home/junjiez1/projects/classes/COMP790/hippo/*/*.vtk'
dataset_fns_list = np.array(glob.glob(dataset_path))

mesh = pv.read(dataset_fns_list[10])
vertices = np.array(mesh.points).astype(np.float32)
faces = mesh.faces.reshape(-1, 4)

mean_vertices = np.zeros(vertices.shape, dtype=np.float32)
for pathname in glob.glob(dataset_path):
    mesh = pv.read(pathname)
    vertices = np.array(mesh.points).astype(np.float32)
    mean_vertices = mean_vertices + vertices

num_subjs = float(len(dataset_fns_list))
mean_vertices = mean_vertices / num_subjs

mean_obj_mesh = pv.PolyData(mean_vertices, faces)
mean_obj_mesh.save('my_mean_obj.vtk')

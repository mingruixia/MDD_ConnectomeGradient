from brainsmash.workbench.geo import volume
from brainsmash.mapgen.sampled import Sampled
import scipy.io as sio

coord_file = "voxel_coordinates.txt"
output_dir = "D:\temp"
filenames = volume(coord_file, output_dir)
brain_map = "g1_Z_T2.txt"
gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'])
surrogate_g1 = gen(n=10000)
sio.savemat('surrogate_g1.mat',{'surrogate_maps':surrogate_g1})
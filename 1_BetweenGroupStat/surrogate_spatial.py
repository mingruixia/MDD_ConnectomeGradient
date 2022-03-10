from brainsmash.workbench.geo import volume
from brainsmash.mapgen.sampled import Sampled
import scipy.io as sio
coord_file = "voxel_coordinates.txt"
output_dir = "D:\Data\SurrogateMapForSpatialCorr"

filenames = volume(coord_file, output_dir)

brain_map = "brain_map_hc_g1.txt"
gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'])
surrogate_maps = gen(n=10000)
sio.savemat('surrogate_maps_hc_g1.mat',{'surrogate_maps':surrogate_maps})


brain_map = "brain_map_hc_g2.txt"
gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'])
surrogate_maps = gen(n=10000)
sio.savemat('surrogate_maps_hc_g2.mat',{'surrogate_maps':surrogate_maps})

brain_map = "brain_map_hc_g3.txt"
gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'])
surrogate_maps = gen(n=10000)
sio.savemat('surrogate_maps_hc_g3.mat',{'surrogate_maps':surrogate_maps})
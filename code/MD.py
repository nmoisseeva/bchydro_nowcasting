
from da_config import *
import numpy as np
from scipy.spatial import KDTree
from util_fcns import plot_domain
import pickle
import matplotlib.pyplot as plt
import os.path
from operator import itemgetter

#==================================Set Up MD class=================================

#define a MD tile class to use in each iteration 
class MDtile:
	def __init__(self,demPoints,demWeights,demDist,mother_idx,HA_idx):
		self.tree_id = mother_idx
		self.coord = demPoints[mother_idx]
		self.S = demWeights[mother_idx]
		self.dist = demDist[mother_idx]
		self.HA = HA_idx
	
	def get_tile(self,demTree,max_dist,landmask_flat):
		tile_dist, tile_id = demTree.query(self.coord,k=9, distance_upper_bound=max_dist)
		tile_dist = tile_dist[1:]							#remove first element (mother point itself)
		trimmed_tile = tile_id[np.isfinite(tile_dist)]		#remove points outside of boundaries
		# if land_flag:
		# 	trimmed_tile = trimmed_tile[landmask_flat[trimmed_tile]==1]						#remove water points
		trimmed_idx = [np.where(tile_id==item)[0][0] for item in trimmed_tile]			#get indecies of used tile points
		trimmed_dist = tile_dist[trimmed_idx]											#get distance to used tile points
		self.tile = trimmed_tile
		self.tile_dist = trimmed_dist


	
	def get_weights(self,dem_flat,params):
		self.tile_weights = np.empty(len(self.tile))
		self.Z = dem_flat[self.tree_id]
		self.haZ = dem_flat[self.HA]
		a,b, Z_ref1, Z_ref2 = params[:] 
		for nCell, cell in enumerate(self.tile):
			cellZ = dem_flat[cell]
			sharing_factor = self.S * ( 1 - ( abs(self.Z - cellZ)/Z_ref1 )**a ) * ( 1 - ( abs(self.haZ - cellZ)/Z_ref2 )**b )
			self.tile_weights[nCell] = sharing_factor

	def update_grid(self,demWeights,demDist,dist_cutoff):
		changed_cells = np.where((self.tile_weights - demWeights[self.tile] > 0.01) & (demDist[self.tile] + self.dist < dist_cutoff) )[0] #update only if the sharing factor is higher and within distance range
		update_ids = self.tile[changed_cells]					#update only if the sharing factor is higher
		update_weights = self.tile_weights[changed_cells]		#update changed tile weights
		update_dist = self.tile_dist[changed_cells]				#update changed distsances
		demWeights[update_ids] = update_weights					#modify tile weights in full array
		demDist[update_ids] = self.dist + update_dist			#modify tile distances in full array
		return demWeights, demDist, update_ids


#==================================Begin Calculations=================================
npy_dir = './npy/'

file_tag = params + (dist_cutoff,)

weights_filename = 'md_weights_%s-%s-%s-%s_r%s' %file_tag
dist_filename = 'md_dist_%s-%s-%s-%s_r%s' %file_tag
id_filename = 'md_id_%s-%s-%s-%s_r%s' %file_tag

basemap_path = basemap_dir + 'basemap_-128.8_-121.0_48.2_51.0_f.pickle'
bm = pickle.load(open(basemap_path,'rb'))   # load the above pickle

md_args = np.load('md_args.npy').item()
elevation, lon_grid, lat_grid, obs, landmask_flat, missing_ids= itemgetter('elevation', 'lon_grid', 'lat_grid', 'obs', 'landmask_flat', 'missing_ids')(md_args)


for key in obs:
	obs[key] = obs[key][missing_ids]

test_lcn = zip(obs['x'],obs['y'])
num_stn = len(test_lcn)
demPoints = zip(lon_grid.ravel(),lat_grid.ravel())
demTree = KDTree(demPoints)
dem_dist, dem_id = demTree.query(test_lcn)
dem_flat = elevation.ravel()

print('Initializing Mother-Daughter routine')
grid_step = np.mean([abs(lon_grid[1,0]-lon_grid[0,0]), abs(lat_grid[0,1]-lat_grid[0,0])])		#horizontal step
max_dist = np.sqrt(2*(grid_step**2))							#max distance for nearst neighbour search (distance to corner grid point)


md_weights_allstn = np.empty((num_stn,len(demPoints)))			#storage arrays for all stns
md_dist_allstn = np.empty((num_stn,len(demPoints)))

for nStn in range(num_stn):		
	print('......... processing station %s of %s' % ((nStn+1),num_stn))						
	demWeights = np.zeros(len(demPoints))						#storage arrays for single stn				
	demDist = np.zeros(len(demPoints))
	HA_idx = dem_id[nStn]										#get index of Honorary Ancestor
	demDist[HA_idx] = 0.000000001 								#set HA distance to ~0 (not actual 0, to avoid later overwrite)
	demWeights[HA_idx] = 1. 									#set HA weight to 1
	mother_idx = dem_id[nStn]									#set mother index (same as HA index for first iteration)
	first_tile = MDtile(demPoints, demWeights, demDist, mother_idx, HA_idx)						#initate an MD tile around the mother
	first_tile.get_tile(demTree,max_dist,landmask_flat) 		#get tile points (which require adjustment)		
	first_tile.get_weights(dem_flat,params)					#get MD weights for the tile
	demWeights, demDist, update_ids = first_tile.update_grid(demWeights, demDist,dist_cutoff) 		#get weights and distances into full array for single stn

	while len(update_ids) > 0:									#as long as there are points to update, continue loop
		modified_list = []										#storage array for modified points
		for nDaughter, daughter in enumerate(update_ids): 		#iterate through modified points
			mother_idx = daughter 								#each point becomes a new mother
			daughter_tile = MDtile(demPoints, demWeights, demDist, mother_idx, HA_idx) 			#repeat MD as above
			daughter_tile.get_tile(demTree,max_dist,landmask_flat)
			daughter_tile.get_weights(dem_flat, params)
			demWeights, demDist, modified_cells = daughter_tile.update_grid(demWeights, demDist, dist_cutoff)
			modified_list.extend(modified_cells) 				#update list of modified points
		modified_list = list(set(modified_list))				#exclude points which were previously calculated from the list
		update_ids = modified_list     							#refresh updated points list

	demWeights[demWeights==0] = np.nan 							#mask uneffected grids
	md_weights_allstn[nStn,:] = demWeights 						#add stn data to full storage array
	demDist[demDist==0] = np.nan 								#repeat for distance
	md_dist_allstn[nStn,:] = demDist

# PYTHON 2.7 BUG: FOR FUTURE IMPLEMENTATION WITH 3.4 ONLY
# md_dict = {'weights':md_weights_allstn, 'dist':md_dist_allstn, 'id':obs['id'][:]}	
# #if a weights file already exists for current configuration append the dictionary with new stations
# if os.path.isfile(weights_path):
# 	print('......... Appending existing weights file')
# 	md_dict_old = np.load(weights_path).item()
# 	for key in md_dict_old.keys():
# 		md_dict[key] = numpy.concatenate((md_dict[key],md_dict_old[key]))
# np.save(weights_path, md_dict)


#CURRENT SAVING SCHEME:
#if a weights file already exists for current configuration append the dictionary with new stations
weights_path = npy_dir + weights_filename +'.npy'
id_path = npy_dir + id_filename +'.npy'
dist_path = npy_dir + dist_filename + '.npy'
if os.path.isfile(weights_path) and os.path.isfile(id_path) and os.path.isfile(dist_path):
	print '......... Appending existing weights file'
	md_weights_old = np.load(weights_path)
	md_id_old = np.load(id_path)
	md_dist_old = np.load(dist_path)
	md_weights_allstn = np.concatenate((md_weights_old,md_weights_allstn))
	md_id = np.concatenate((md_id_old,obs['id'][:]))
	md_dist_allstn = np.concatenate((md_dist_old,md_dist_allstn))
else:
	md_id = obs['id'][:]
np.save(weights_path, md_weights_allstn)
np.save(dist_path, md_dist_allstn)
np.save(id_path, md_id)

print('......... MD weights saved in: ' + weights_path)




#==================================Plotting=================================
#plot the weights over topography
bm_fig = plot_domain(bm)
bm.imshow(elevation.T, cmap=plt.cm.gist_yarg, alpha = 0.5, origin='upper' )
bm.drawcoastlines()
plt.title('Station Locations and MD Weights', fontweight='bold')
for nStn in range(num_stn):
	temp = np.reshape(md_weights_allstn[nStn,:], np.shape(elevation))
	bm.imshow(temp.T, origin='upper', alpha = 0.5)
plt.colorbar(label='sharing factor',orientation='horizontal')
plt.savefig(fig_dir + weights_filename +'.pdf')
plt.close()

#plot distances over topography
bm_fig = plot_domain(bm)
bm.imshow(elevation.T, cmap=plt.cm.gist_yarg, alpha = 0.5, origin='upper' )
plt.title('Station Locations and MD Anisotropic Distance', fontweight='bold')
for nStn in range(num_stn):
	temp = np.reshape(md_dist_allstn[nStn,:], np.shape(elevation))
	bm.imshow(temp.T, origin='upper', alpha = 0.5)
plt.colorbar(label='anisotropic distance [deg]',orientation='horizontal')
bm.scatter(np.array(obs['x']),np.array(obs['y']),s=6,marker='o', color ='r')
plt.savefig(fig_dir + 'MD_dist_%s-%s-%s-%s_r%s.pdf' %file_tag)
plt.close()









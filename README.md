# Overview

Fractal Dimension Neighbor Embedding (fdNE) algorithm in application to fitness of deletions in virus genome across groups (passages). 

The genome deletions (junctions) are defined by two coordinates (Start, End), frequency (psi1), standard deviation of the frequency (sd1) and can be depicted as points on the plane—genome positions on both coordinate axes.

The main script named fdNE.py. 
It performs two steps: 
1. Finding neighborhoods of deletions;
2. Calculating of significance of each neighborhoods throught all groups (passages).

## How to use 
```
Optional arguments fdNE.py:
  -h, --help            show this help message and exit
  -d DEL_TABLES, --del_tables DEL_TABLES
                       Comma seaprated deleteions tables paths. File names must contain sample number and group (passage) number. e.g. "vero_sample1_group2.txt,vero_sample1_group3.txt"
  -f FREQ_FILT, --freq_filt FREQ_FILT
                        First threshold for frequency filtering. Points with frequency more than threshold will be used for neighborhood construction. (default=auto)
  -s FREQ_FILT_2, --freq_filt_2 FREQ_FILT_2
                        Second threshold for frequency filtering. Points with frequency more than threshold will be added to neighborhood analysis after their construction. (default=0.0)
  -l DEL_LENGTH, --del_length DEL_LENGTH
                        Deletions length threshold. (default=400)
  -c, --only_coords     Coords print flag.
  --do_two_step         Do two step neighborhood identification?.
  --do_mod_step         Do modified step neighborhood identification?.
  
  -K NUM_K_POINTS, --num_K_points NUM_K_POINTS
                        Number of K random points. (default=5)
  -A A_VALUE, --A_value A_VALUE
                        Value for rounding distances between all points and K
                        points. (default=auto)
```
### Output files description
Output will write to ./out_deletions_table.
The number of files may vary because the number of neighborhoods may vary.
```
all_points_frq_<freq_thresh_1>.txt - coordinates of point with frequencies (psi1) more than <freq_thresh_1>.
all_points_frq_<freq_thresh_2>.txt - coordinates of point with frequencies (psi1) more than <freq_thresh_2>.
1590.0_2810.0_coords_points - coordinates of points (with frequencies (psi1) more than  <freq_thresh_1> ) which belong to neighborhood  with center coord 1590.0 2810.0
cov_junct_per_pos_out_neighb_info.txt - deletions average coverage in every group of samples. For all neighborhoods.
num_points_out_neighb_info.txt - number of points in every group of samples for all neighborhoods.
sum_significance_out_neighb_info.txt - sum of differences of each SignPsi1  minus AvSignJunct in every group of samples for all neighborhoods. singPsi1 = psi1/sd1 calculated for each point, AvSignJunct - average  singPsi1 for all points of group.
significance_out_neighb_info.txt - siginificance in every group of samples for all neighborhoods.  
constructed_neighborhoods_new_centers_step_2.txt - file with information about costructed  neighborhoods included points with frequencies (psi1) more than <freq_thresh_1> (center of neighborhood, distances from all pointsof  neighborhood to center, pvalues, dimensions, all point which belong to neighborhood).
constructed_clusters_step_2_with_all_points.txt - file with information about costructed clusters without inersections included points with frequencies (psi1) more than <freq_thresh_2>
```

### Example
Example input files in /example directory.

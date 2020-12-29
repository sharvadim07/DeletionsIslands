import math
import random
import numpy
from scipy.spatial import distance
from abc import ABC, abstractmethod
from scipy.stats import poisson

import neighborhood

class Point(ABC):
    @abstractmethod
    def get_value(self):
        pass

    @abstractmethod
    def calc_dist(self, Point):
        pass

class DistanceMatrix():
    def __init__(self, points_tuple):
        #self.dist_matrix = [[0.0 for x in range(len(points_tuple))] for y in range(len(points_tuple))]
        self.size = len(points_tuple)
        import multiDim
        if type(points_tuple[0]) is multiDim.MultiDimPoint:
            param_vectors = numpy.asarray(
                [point.param_tuple for point in points_tuple])
            self.dist_matrix = distance.cdist(
                param_vectors, param_vectors, 'euclidean')
            numpy.fill_diagonal(self.dist_matrix, 1E-299)
        else:
            self.dist_matrix = numpy.zeros(
                (len(points_tuple), len(points_tuple)), float)
            for i in range(0, len(points_tuple)):
                for j in range(i, len(points_tuple)):
                    dist = points_tuple[i].calc_dist(points_tuple[j])
                    # if dist == 0.0:
                    #     dist = 2E-19
                    self.dist_matrix[i, j] = dist
                    self.dist_matrix[j, i] = dist

def add_all_points(points_tuple, points_tuple_no_filter, cur_space):
    # Add other points without filtering by frequency and length
    for neighb in cur_space.neighb_sort_by_feature_list:
        neighb.add_points_no_filter_to_neighb(
            points_tuple, points_tuple_no_filter)
    # Resorting neighborhoods in space
    cur_space.sort_neighbs_by_feature()
    new_points_tuple = tuple(points_tuple + points_tuple_no_filter)

    cur_space.points_tuple = new_points_tuple
    # TODO: Can off dist_matrix recalc when need speed of calculations
    #cur_space.dist_matrix = space.DistanceMatrix(cur_space.points_tuple)
    cur_space.dist_matrix = None
    return new_points_tuple

def get_points_tuple_in_cur_neghbs(cur_space):
    all_neighbs_closest_points = [neighb.closest_points + [[neighb.center_point_ind, 0.0]]
                                  for neighb in cur_space.neighb_sort_by_feature_list]
    all_neighbs_points_ind_list = []
    for closest_points in all_neighbs_closest_points:
        for point_ind_dist in closest_points:
            all_neighbs_points_ind_list.append(point_ind_dist)
        points_tuple_in_cur_neighbs = tuple(dict.fromkeys(
            [cur_space.points_tuple[point_ind_dist[0]]
             for point_ind_dist in all_neighbs_points_ind_list]))
    return points_tuple_in_cur_neighbs

def transform_dist_matrix(dist_matrix, neighborhood_list):
    import copy
    trans_dist_matrix = copy.deepcopy(dist_matrix)
    for neighb in neighborhood_list:
        pvalue_str = ''
        if neighb.pvalue != 0.0:
            pvalue_str = '%.2E' % neighb.pvalue
        else:
            pvalue_str = '%.2E' % float('1.0E-300')
        number = float(pvalue_str.split('E-')[0])
        power = int(pvalue_str.split('E-')[1])
        factor = number*10.0**(-math.log(power))
        for point_1_ind_dist in neighb.closest_points:
            for point_2_ind_dist in neighb.closest_points:
                trans_dist_matrix.dist_matrix[point_1_ind_dist[0], point_2_ind_dist[0]
                                              ] = dist_matrix.dist_matrix[point_1_ind_dist[0], point_2_ind_dist[0]]*factor
    return trans_dist_matrix

def print_dist_mat(dist_matrix, out_file_name):
    import os
    with open(os.path.basename(out_file_name), 'w') as out_file:
        for i in range(dist_matrix.size):
            line_matrix = ['%.5f' %
                           dist for dist in dist_matrix.dist_matrix[i]]
            out_file.write('\t'.join(line_matrix) + '\n')

class Space():
    def __init__(self, points_tuple,
                 num_K_points,
                 A,
                 zone_size_percent=50,
                 Q=1000,
                 small_circle_size=7,
                 neighb_post_process=True,
                 second_step=False,
                 num_of_random_centers_from_isalnd=1,
                 sort_neighb_by_density=False,
                 pvalue_threshold=0.0000001,
                 min_neighb_part_from_all=0.01):
        self.points_tuple = points_tuple
        self.dist_matrix = DistanceMatrix(points_tuple)
        if A == 'auto':
            self.A = self.dist_matrix.dist_matrix.max()/100
            if self.A < 1:
                self.A = 1
        else:
            self.A = float(A)
        self.num_K_points = num_K_points
        self.min_neighb_size = small_circle_size
        self.sort_neighb_by_density = sort_neighb_by_density
        self.K_points_tuple = self.select_K_points()
        self.pvalue_points_in_neighb_dict = None
        #self.hash_islands_dict = self.hash_islands(A, min_island_size=7)
        if min_neighb_part_from_all*len(points_tuple) > self.min_neighb_size:
            self.min_neighb_size = min_neighb_part_from_all*len(points_tuple)
        self.hash_islands_dict = self.hash_islands(
            self.A, min_island_size=self.min_neighb_size)
        self.neighborhood_list = self.create_neighborhoods(
            num_of_random_centers_from_isalnd)
        # Filter by size and sort
        self.neighb_sort_by_feature_list = \
            [neighb for neighb in self.detect_most_enriched_neighborhood(zone_size_percent, Q, second_step)
                if neighb.size >= self.min_neighb_size]
        self.sort_neighbs_by_feature()
        self.neighborhood_list = []

        # Debug
        # import os
        # cwd = os.getcwd()
        # os.makedirs('out_debug', exist_ok=True)
        # os.chdir('out_debug')
        # for neighb in self.neighb_sort_by_feature_list:
        #     with open(str(neighb.center_point_ind) + '_test_debug.txt', 'w') as out_file:
        #         out_file.write('\t'.join(list(map(str, self.points_tuple[neighb.center_point_ind].param_tuple))) + '\n')
        #         for i, pval in enumerate(neighb.pvalue_list):
        #             pval_thresh = 1/(i+4)
        #             out_file.write(str(i) + '\t' + str(pval) + '\t' + str(pval_thresh) + '\n')
        # exit()

        if len(self.neighb_sort_by_feature_list) == 0:
            raise RuntimeError(
                'Neighborhoods not found with current input data!')

        # self.find_new_centers_of_neighborhoods()

        if neighb_post_process:
            # Find best pvalues for point located in found neighborhoods
            self.find_pvalue_for_points_in_all_neighb(small_circle_size)
            self.pvalue_points_in_neighb_med = \
                numpy.median(list(self.pvalue_points_in_neighb_dict.values()))
            # Delete points with low lambda value from all neighborhood
            # self.del_low_lambda_points_from_all_neighb()

            # Resolving intersection in all cases (multi dim and two dim)
            # Find intersection and delete points from neighborhood
            # with smaller pvalues or densities
            if sort_neighb_by_density:
                self.get_neighborhoods_without_intersection_by_density()
            else:
                self.get_neighborhoods_without_intersection_by_pval()
            # import multiDim
            # if type(self.points_tuple[0]) is multiDim.MultiDimPoint:
            #     self.get_neighborhoods_without_intersection(Q)

            # Filter points by pvalue
            self.filter_points_by_pvalue(self.pvalue_points_in_neighb_med)

    def set_points_tuple(self, points_tuple):
        self.points_tuple = points_tuple
        self.dist_matrix = DistanceMatrix(points_tuple)

    def find_new_centers_of_neighborhoods(self):
        for neighb in self.neighb_sort_by_feature_list:
            neighb.find_new_center(self.dist_matrix)

    def sort_neighbs_by_feature(self):
        if self.sort_neighb_by_density:
            self.neighb_sort_by_feature_list = \
                sorted(self.neighb_sort_by_feature_list,
                       key=lambda x: x.density, reverse=True)
        else:
            self.neighb_sort_by_feature_list = \
                sorted(self.neighb_sort_by_feature_list,
                       key=lambda x: x.pvalue)

    def select_K_points(self):
        i = self.num_K_points
        K_points_list = []
        while i != 0:
            r_num = random.randint(0, len(self.points_tuple) - 1)
            if r_num not in K_points_list:
                if len(K_points_list) != 0:
                    dist_to_prev = self.dist_matrix.dist_matrix[K_points_list[-1], r_num]
                    if dist_to_prev < numpy.max(self.dist_matrix.dist_matrix)/10:
                        continue
                K_points_list.append(r_num)
                i -= 1
        return tuple(K_points_list)

    def rounded_dist_to_K(self, A=10):
        #rounded_dist_to_K_mat = [[1.0 for x in range(len(self.K_points_tuple))] for y in range(len(self.points_tuple))]
        rounded_dist_to_K_mat = numpy.full(
            (len(self.points_tuple), len(self.K_points_tuple)), 1.0, float)
        for i in range(len(self.points_tuple)):
            for j, K in enumerate(self.K_points_tuple):
                if i != K:
                    rounded_dist_to_K_mat[i, j] = round(
                        self.dist_matrix.dist_matrix[i, K]/A, 0)
                    # if rounded_dist_to_K_mat[i,j] == 0.0:
                    #     rounded_dist_to_K_mat[i,j] = 1.0
        return rounded_dist_to_K_mat

    def my_hash(self, x_ind, rd):
        hash_res = 1.0
        for i in range(len(self.K_points_tuple)):
            # hash_el = math.log10(rd[x_ind][i])
            hash_el = math.log10(rd[x_ind, i] + 1.0)
            if hash_el == 0.0:
                hash_el = 1.0
            hash_res *= hash_el
        return round(hash_res, 2)

    def hash_islands(self, A, min_island_size=7):
        hash_dict = {}
        rd = self.rounded_dist_to_K(A)
        for i in range(len(self.points_tuple)):
            # TODO: Check
            #hash_res = str(self.my_hash(i, rd)).split('e')[1]
            hash_res = self.my_hash(i, rd)
            if hash_res in hash_dict:
                hash_dict[hash_res].append(i)
            else:
                hash_dict[hash_res] = [i]
        islands_to_remove_list = []
        for island in hash_dict:
            if len(hash_dict[island]) < min_island_size:
                islands_to_remove_list.append(island)
        for island in islands_to_remove_list:
            hash_dict.pop(island)
        return hash_dict

    def create_neighborhoods(self, num_of_random_centers_from_isalnd=1):
        neighborhood_list = []
        for island in self.hash_islands_dict:
            for i in range(num_of_random_centers_from_isalnd):
                neighb_center_ind = self.hash_islands_dict[island][random.randint(
                    0, len(self.hash_islands_dict[island]) - 1)]
                import twoDim
                import multiDim
                if type(self.points_tuple[0]) is twoDim.TwoDimPointDel:
                    neighb = twoDim.NeighborhoodTwoDim(neighb_center_ind)
                elif type(self.points_tuple[0]) is multiDim.MultiDimPoint:
                    neighb = multiDim.NeighborhoodMultiDim(neighb_center_ind)
                else:
                    raise TypeError("Unknown type in points tuple in space!")
                neighb.init_neighborhood(
                    dist_matrix=self.dist_matrix, init_size=3)
                neighborhood_list.append(neighb)
        return neighborhood_list

    def detect_most_enriched_neighborhood(self, zone_size_percent, Q, second_step=False):
        for neighb in self.neighborhood_list:
            import twoDim
            import multiDim
            if type(self.points_tuple[0]) is twoDim.TwoDimPointDel:
                # neighb.find_enrichment_area_radius_2_steps(Q)
                if not second_step:
                    neighb.find_enrichment_area_radius(zone_size_percent, Q)
                else:
                    neighb.find_enrichment_area_radius_2_steps(
                        zone_size_percent, Q)
                neighb.dimension = 2
                neighb.volume = math.pi*(neighb.r[1]**2)
                neighb.density = neighb.size/neighb.volume
            elif type(self.points_tuple[0]) is multiDim.MultiDimPoint:
                if not second_step:
                    neighb.find_enrichment_area_radius(zone_size_percent, Q)
                else:
                    neighb.find_enrichment_area_radius_2_steps(
                        zone_size_percent, Q)
            else:
                raise TypeError("Unknown type in points tuple in space!")
        return sorted(self.neighborhood_list, key=lambda x: x.pvalue)

    def del_point_from_neighb(self,
                              neighb_ind_from,
                              neighb_ind_to,
                              points_intersects,
                              neighborhood_to_del_list,
                              neighb_intersections_dict,
                              full_intersection_flag=False):
        # Check when full intersection flag always False
        #full_intersection_flag = False
        if self.neighb_sort_by_feature_list[neighb_ind_from].size - len(points_intersects) < self.min_neighb_size:
            neighborhood_to_del_list.append(neighb_ind_from)
        else:
            neighb_ind_from_size = self.neighb_sort_by_feature_list[neighb_ind_from].size
            # Simple points removing
            for point_to_remove in points_intersects:
                self.neighb_sort_by_feature_list[neighb_ind_from].remove_point(
                    point_to_remove)
            # Saving deletion history of points followed by recalculation of the PValue
            if neighb_ind_from not in neighb_intersections_dict:
                # Save initial size
                neighb_intersections_dict[neighb_ind_from] = [
                    neighb_ind_from_size]
                neighb_intersections_dict[neighb_ind_from].append(
                    (neighb_ind_to, len(points_intersects), full_intersection_flag)
                )
            else:
                neighb_intersections_dict[neighb_ind_from].append(
                    (neighb_ind_to, len(points_intersects), full_intersection_flag)
                )

    def neighb_pvalue_recalculation(self, neighb_intersections_dict, neighborhood_to_del_list):
        if 0 not in neighb_intersections_dict:
            self.neighb_sort_by_feature_list[0].pvalue = neighborhood.pval_calc(self.neighb_sort_by_feature_list[0].lambda_val,
                                                                                self.neighb_sort_by_feature_list[0].volume,
                                                                                self.neighb_sort_by_feature_list[0].size)
        for neighb_from_ind in range(0, len(self.neighb_sort_by_feature_list)):
            if neighb_from_ind not in neighborhood_to_del_list and neighb_from_ind in neighb_intersections_dict:
                neighb_from_size = neighb_intersections_dict[neighb_from_ind][0]
                neighb_from_lambda = self.neighb_sort_by_feature_list[neighb_from_ind].lambda_val
                neighb_from_volume = \
                    self.neighb_sort_by_feature_list[neighb_from_ind].volume
                # For torus
                # Find Volume of torus
                deleted_torus_size = \
                    sum([info[1] for info in neighb_intersections_dict[neighb_from_ind][1:]
                         if (info[2] == True and info[0] not in neighborhood_to_del_list)])
                Vtorus = neighb_from_volume
                new_neighb_from_size = neighb_from_size
                if deleted_torus_size != 0:
                    sum_of_deleted_volumes = \
                        sum([self.neighb_sort_by_feature_list[info[0]].r[1]**self.neighb_sort_by_feature_list[neighb_from_ind].dimension
                             for info in neighb_intersections_dict[neighb_from_ind][1:]
                             if (info[2] == True and info[0] not in neighborhood_to_del_list)])
                    Vtorus = neighb_from_volume - sum_of_deleted_volumes
                    # Update pvalue after torus deletion
                    # self.neighb_sort_by_feature_list[neighb_from_ind].pvalue = \
                    #     (neighb_from_lambda*Vtorus)**(neighb_from_size - deleted_torus_size)*math.e**(-neighb_from_lambda*Vtorus)/ \
                    #         math.factorial(neighb_from_size - deleted_torus_size)
                    self.neighb_sort_by_feature_list[neighb_from_ind].volume = Vtorus
                    self.neighb_sort_by_feature_list[neighb_from_ind].density = \
                        self.neighb_sort_by_feature_list[neighb_from_ind].size/Vtorus
                    self.neighb_sort_by_feature_list[neighb_from_ind].pvalue = neighborhood.pval_calc(neighb_from_lambda,
                                                                                                      Vtorus,
                                                                                                      neighb_from_size - deleted_torus_size)
                    new_neighb_from_size = neighb_from_size - deleted_torus_size
                # For sickle
                deleted_sickle_size = sum([info[1] for info in neighb_intersections_dict[neighb_from_ind][1:]
                                           if (info[2] != True and info[0] not in neighborhood_to_del_list)])
                if deleted_sickle_size != 0:
                    Vsickle = (new_neighb_from_size -
                               deleted_sickle_size)*Vtorus/new_neighb_from_size
                    # Update pvalue after sickle deletion
                    # self.neighb_sort_by_feature_list[neighb_from_ind].pvalue = \
                    #     (neighb_from_lambda*Vsickle)**(new_neighb_from_size - deleted_sickle_size)*math.e**(-neighb_from_lambda*Vsickle)/ \
                    #         math.factorial(new_neighb_from_size - deleted_sickle_size)
                    self.neighb_sort_by_feature_list[neighb_from_ind].volume = Vsickle
                    self.neighb_sort_by_feature_list[neighb_from_ind].density = \
                        self.neighb_sort_by_feature_list[neighb_from_ind].size/Vsickle
                    self.neighb_sort_by_feature_list[neighb_from_ind].pvalue = neighborhood.pval_calc(neighb_from_lambda,
                                                                                                      Vsickle,
                                                                                                      neighb_from_size - deleted_sickle_size)
                    if (new_neighb_from_size - deleted_sickle_size) != self.neighb_sort_by_feature_list[neighb_from_ind].size:
                        raise ValueError(
                            "Neighb size after deletion of points and after pvalue recalculation not match!")

    def get_neighborhoods_without_intersection_by_density(self, only_full_intersection = False):
        # Dict with information about
        # deferred deletion of points from neighb
        # keys - indexes of neighb in neighb_sort_by_feature_list
        neighb_intersections_dict = {}
        # Find neighborhood which not have intersection
        neighborhood_to_del_list = []
        for i, neighb in enumerate(self.neighb_sort_by_feature_list[:-1]):
            if i in neighborhood_to_del_list:
                continue
            neighb_closest_point_ind_tuple = tuple(closest_point_ind_dist[0]
                                                   for closest_point_ind_dist in
                                                   neighb.closest_points + [[neighb.center_point_ind, 0]])
            for j, neighb_next in enumerate(self.neighb_sort_by_feature_list[i+1:]):
                if (i + j + 1) in neighborhood_to_del_list:
                    continue
                neighb_next_points_intersects = []
                # Loop for find Points which intersects neighb and neighb_next
                for neighb_next_closest_point in neighb_next.closest_points + [[neighb_next.center_point_ind, 0]]:
                    if neighb_next_closest_point[0] in neighb_closest_point_ind_tuple:
                        neighb_next_points_intersects.append(
                            neighb_next_closest_point[0])
                if len(neighb_next_points_intersects) != 0:
                    if (neighb_next.size - len(neighb_next_points_intersects)) < self.min_neighb_size and \
                            i + j + 1 not in neighb_intersections_dict:
                        neighborhood_to_del_list.append(i + j + 1)
                    elif only_full_intersection == False:
                        # TODO: Maybe need do same as
                        # Implement intersection cut smaller neighb
                        # from bigger neighb with similar pvalue
                        if (neighb.size - len(neighb_next_points_intersects)) < self.min_neighb_size  and \
                                i not in neighb_intersections_dict:
                            self.del_point_from_neighb(i + j + 1,
                                                       i,
                                                       neighb_next_points_intersects,
                                                       neighborhood_to_del_list,
                                                       neighb_intersections_dict,
                                                       True)
                        else:
                            self.del_point_from_neighb(i + j + 1,
                                                       i,
                                                       neighb_next_points_intersects,
                                                       neighborhood_to_del_list,
                                                       neighb_intersections_dict)
        # PValue recalculation
        # keep in mind that some areas have been removed
        if only_full_intersection == False:
            self.neighb_pvalue_recalculation(
                neighb_intersections_dict, neighborhood_to_del_list)
        for neighborhood_to_del in sorted(neighborhood_to_del_list, reverse=True):
            del self.neighb_sort_by_feature_list[neighborhood_to_del]
        # Reorder by new density
        self.neighb_sort_by_feature_list = sorted(
            self.neighb_sort_by_feature_list, key=lambda x: x.density, reverse=True)
        # TODO: Can be turn off
        # self.find_new_centers_of_neighborhoods()

    def get_neighborhoods_without_intersection_by_pval(self, pval_similarity=0.8):
        # Dict with information about
        # deferred deletion of points from neighb
        # keys - indexes of neighb in neighb_sort_by_feature_list
        neighb_intersections_dict = {}
        # Find neighborhood which not have intersection
        neighborhood_to_del_list = []
        for i, neighb in enumerate(self.neighb_sort_by_feature_list[:-1]):
            if i in neighborhood_to_del_list:
                continue
            neighb_closest_point_ind_tuple = tuple(closest_point_ind_dist[0]
                                                   for closest_point_ind_dist in
                                                   neighb.closest_points + [[neighb.center_point_ind, 0]])
            neighb_points_intersects = []
            for j, neighb_next in enumerate(self.neighb_sort_by_feature_list[i+1:]):
                if (i + j + 1) in neighborhood_to_del_list:
                    continue
                neighb_next_points_intersects = []
                # Loop for find Points which intersects neighb and neighb_next
                for neighb_next_closest_point in neighb_next.closest_points + [[neighb_next.center_point_ind, 0]]:
                    if neighb_next_closest_point[0] in neighb_closest_point_ind_tuple:
                        neighb_next_points_intersects.append(
                            neighb_next_closest_point[0])
                if len(neighb_next_points_intersects) != 0:
                    if (neighb_next.size - len(neighb_next_points_intersects)) < self.min_neighb_size and \
                            i + j + 1 not in neighb_intersections_dict:
                        # Implement intersection cut smaller neighb
                        # from bigger neighb with similar pvalue
                        if (neighb_next.size <= neighb.size/2 and neighb.pvalue/neighb_next.pvalue >= pval_similarity):
                            neighb_points_intersects = neighb_next_points_intersects
                            neighb_next_points_intersects = []
                            # Delete points from neighb with better pvalue
                            self.del_point_from_neighb(i,
                                                       i + j + 1,
                                                       neighb_points_intersects,
                                                       neighborhood_to_del_list,
                                                       neighb_intersections_dict,
                                                       True)
                        else:
                            neighborhood_to_del_list.append(i + j + 1)
                    else:
                        # TODO: Maybe need do same as
                        # Implement intersection cut smaller neighb
                        # from bigger neighb with similar pvalue
                        if (neighb.size - len(neighb_next_points_intersects)) < self.min_neighb_size and \
                                i not in neighb_intersections_dict:
                            self.del_point_from_neighb(i + j + 1,
                                                       i,
                                                       neighb_next_points_intersects,
                                                       neighborhood_to_del_list,
                                                       neighb_intersections_dict,
                                                       True)
                        else:
                            self.del_point_from_neighb(i + j + 1,
                                                       i,
                                                       neighb_next_points_intersects,
                                                       neighborhood_to_del_list,
                                                       neighb_intersections_dict)
        # PValue recalculation
        # keep in mind that some areas have been removed
        self.neighb_pvalue_recalculation(
            neighb_intersections_dict, neighborhood_to_del_list)
        for neighborhood_to_del in sorted(neighborhood_to_del_list, reverse=True):
            del self.neighb_sort_by_feature_list[neighborhood_to_del]
        # Reorder by new pvalue
        self.neighb_sort_by_feature_list = sorted(
            self.neighb_sort_by_feature_list, key=lambda x: x.pvalue)
        # TODO: Can be turn off
        # self.find_new_centers_of_neighborhoods()

    def del_bad_pvalue_points_from_all_neighb(self):
        for neighb in self.neighb_sort_by_feature_list:
            neighb.del_points_with_bad_pvalue(self.dist_matrix)

    def find_pvalue_for_points_in_all_neighb(self, small_circle_size=7):
        self.pvalue_points_in_neighb_dict = {}
        self.small_circles_dict = {}
        for neighb in self.neighb_sort_by_feature_list:
            # If neighb size equal or less than small_circle_size
            # then set pvalues for all points in this neighb
            # same as neighb pvalue
            # if neighb.size <= small_circle_size:
            #    pvalue = neighb.pvalue
            #    for point_i_dist in neighb.closest_points + [[neighb.center_point_ind, 0.0]]:
            #        if point_i_dist[0] not in self.pvalue_points_in_neighb_dict or \
            #            pvalue < self.pvalue_points_in_neighb_dict[point_i_dist[0]]:
            #            self.pvalue_points_in_neighb_dict[point_i_dist[0]] = pvalue
            if neighb.size < small_circle_size:
                raise ValueError(
                    'Neighborhood size less than small circle size!')
            else:
                for point_i_dist in neighb.closest_points + [[neighb.center_point_ind, 0.0]]:
                    small_circle_points_list = \
                        neighb.find_closest_points_to_point(
                            self.dist_matrix, point_i_dist[0], inside_neighb=True)[:small_circle_size - 1]
                    pvalue = neighb.calc_pval_for_cur_size_small_circle(
                        small_circle_points_list, small_circle_size)
                    if point_i_dist[0] not in self.pvalue_points_in_neighb_dict or \
                            pvalue < self.pvalue_points_in_neighb_dict[point_i_dist[0]]:
                        self.pvalue_points_in_neighb_dict[point_i_dist[0]] = pvalue
                        small_circle = neighborhood.Neighborhood(
                            point_i_dist[0])
                        small_circle.init_neighborhood(closest_points=small_circle_points_list,
                                                       r=small_circle_points_list[-1],
                                                       pvalue=pvalue,
                                                       dimension=neighb.dimension,
                                                       lambda_val=neighb.lambda_val,
                                                       min_r=0.0)
                        small_circle.size = len(small_circle_points_list) + 1
                        self.small_circles_dict[small_circle.center_point_ind] = small_circle
        # return pvalue_points_in_neighb_dict, small_circles_dict

    def filter_points_by_pvalue(self, pvalue_threshold=0.000001):
        neigb_to_del_list = []
        for i, neighb in enumerate(self.neighb_sort_by_feature_list):
            if i in neigb_to_del_list:
                continue
            for point_i_dist in neighb.closest_points + [[neighb.center_point_ind, 0.0]]:
                if self.pvalue_points_in_neighb_dict[point_i_dist[0]] > pvalue_threshold:
                    neighb.remove_point(point_i_dist[0])
                    if point_i_dist[0] in self.small_circles_dict:
                        del self.small_circles_dict[point_i_dist[0]]
                    if neighb.size < self.min_neighb_size and i not in neigb_to_del_list:
                        neigb_to_del_list.append(i)
                    # TODO: Not worked!!!
                    # elif neighb.has_center == False:
                    #     neighb.center_point_ind = neighb.closest_points[0][0]
        for neigb_to_del in sorted(neigb_to_del_list, reverse=True):
            del self.neighb_sort_by_feature_list[neigb_to_del]

    # Small circles graphs searching
    def find_small_circles_graphs(self, points_tuple_no_filter, down_quartile_lambda = None):
        small_circles_graphs_dict = {}
        if down_quartile_lambda == None:
            lambda_val_list = [
                neighb.lambda_val for neighb in self.neighb_sort_by_feature_list]
            # self.all_neighbs_avg_lambda = \
            #     numpy.average(lambda_val_list)
            quartile_idx = int((len(lambda_val_list)-1)/4)
            down_quartile_lambda = \
                sorted(lambda_val_list)[quartile_idx]
        # work_small_circles_dict will decrease    
        work_small_circles_dict = self.small_circles_dict.copy()
        points_ind_set_added = set()
        while len(work_small_circles_dict) != 0:
            del_small_circle_from_work_list = []
            first_small_circle_center_ind = next(iter(work_small_circles_dict))
            small_circles_graphs_dict[first_small_circle_center_ind] = \
                SmallCirclesGraph(
                    work_small_circles_dict[first_small_circle_center_ind])
            for small_circle_cnt_ind in work_small_circles_dict:
                for small_circles_graph in small_circles_graphs_dict:
                    res_add = small_circles_graphs_dict[small_circles_graph].add(self.small_circles_dict[small_circle_cnt_ind],
                                                                                 self.dist_matrix,
                                                                                 down_quartile_lambda,
                                                                                 self.points_tuple,
                                                                                 points_tuple_no_filter,
                                                                                 points_ind_set_added,
                                                                                 len(work_small_circles_dict))
                    if res_add != 2:
                        if res_add != 1:
                            points_ind_set_added = res_add
                        del_small_circle_from_work_list.append(
                            small_circle_cnt_ind)
                        break
            for del_small_circle_from_work in del_small_circle_from_work_list:
                del work_small_circles_dict[del_small_circle_from_work]

        # Present graphs as neighborhoods
        self.neighb_sort_by_feature_list = []
        self.points_tuple = tuple(self.points_tuple + points_tuple_no_filter)
        for small_circles_graph in small_circles_graphs_dict:
            if len(small_circles_graphs_dict[small_circles_graph].points_set) > 0:
                neighb_small_circles_graph = neighborhood.Neighborhood(
                    small_circles_graph)
                neighb_small_circles_graph.init_neighborhood(
                    closest_points=[[point_ind, -1] for point_ind in small_circles_graphs_dict[small_circles_graph].points_set])
                neighb_small_circles_graph.size = len(
                    neighb_small_circles_graph.closest_points) + 1
                self.neighb_sort_by_feature_list.append(
                    neighb_small_circles_graph)
        pass
    
class SmallCirclesGraph():
    def __init__(self, small_circle_init):
        self.small_circles_dict = {
            small_circle_init.center_point_ind: small_circle_init}
        self.small_circles_pairs_dict = {}
        self.points_set = set()

    def add(self,
            new_small_circle,
            dist_matrix,
            base_lambda,
            points_tuple,
            points_tuple_no_filter,
            points_ind_set_added,
            small_circles_dict_size):
        if new_small_circle.center_point_ind in self.small_circles_dict:
            return 1
        # Init new neighb from two small circles
        neighb_from_two_small_circle = neighborhood.Neighborhood(
                new_small_circle.center_point_ind)
        for small_circle_cnt_ind in self.small_circles_dict:
            cur_small_circle = self.small_circles_dict[small_circle_cnt_ind]
            neighb_from_two_small_circle_dimension = -1
            neighb_from_two_small_circle_lambda = base_lambda
            if cur_small_circle.dimension > new_small_circle.dimension:
                neighb_from_two_small_circle_dimension = new_small_circle.dimension
                #neighb_from_two_small_circle_lambda = new_small_circle.lambda_val
            else:
                neighb_from_two_small_circle_dimension = cur_small_circle.dimension
                #neighb_from_two_small_circle_lambda = cur_small_circle.lambda_val

            # Add not duplicated points
            comm_points = [point_ind_dist[0] for point_ind_dist in cur_small_circle.closest_points + new_small_circle.closest_points]
            comm_points.append(cur_small_circle.center_point_ind)
            comm_points.append(new_small_circle.center_point_ind)
            comm_points = set(comm_points)

            neighb_from_two_small_circle_points = [[point_ind, 0.0] for point_ind in comm_points]
            neighb_from_two_small_circle.center_point_ind = neighb_from_two_small_circle_points[-1][0]
            neighb_from_two_small_circle.init_neighborhood(
                closest_points=neighb_from_two_small_circle_points[:-1],
                lambda_val=neighb_from_two_small_circle_lambda,
                dimension=neighb_from_two_small_circle_dimension,
                min_r = 0.0)
            neighb_from_two_small_circle.find_new_center(dist_matrix)

            # Find pvalue for neighb_from_two_small_circle
            pvalue = poisson.sf(neighb_from_two_small_circle.size,
                                neighb_from_two_small_circle.lambda_val *
                                (neighb_from_two_small_circle.r[1]**neighb_from_two_small_circle.dimension))
            #pvalue_threshold = 1/(small_circles_dict_size*(small_circles_dict_size-1)/2)
            pvalue_threshold = 0.05
            if pvalue < pvalue_threshold:
                neighb_from_two_small_circle.pvalue = pvalue

                # Add all points
                neighb_from_two_small_circle.add_points_no_filter_to_neighb(points_tuple,
                                                                            points_tuple_no_filter,
                                                                            points_ind_set_added)
                cur_added_list = [ind_dist[0] for ind_dist in neighb_from_two_small_circle.closest_points] + \
                    [neighb_from_two_small_circle.center_point_ind]
                points_ind_set_added.update(cur_added_list)

                self.points_set.update(cur_added_list)

                self.small_circles_dict[new_small_circle.center_point_ind] = new_small_circle
                self.small_circles_pairs_dict[(small_circle_cnt_ind, new_small_circle.center_point_ind)] = \
                    neighb_from_two_small_circle
                return points_ind_set_added
        return 2

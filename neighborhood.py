import math
from scipy.stats import poisson
import scipy.optimize


class Neighborhood():
    def __init__(self, point_center_ind):
        self.center_point_ind = point_center_ind
        self.has_center = True

    def init_neighborhood(self,
                          init_size=3,
                          closest_points=[],
                          dist_matrix=None,
                          pvalue=1.1,
                          dimension=0.0,
                          r=[-1, 0.0],
                          lambda_val=-1,
                          volume=-1,
                          density=-1,
                          min_r=-1):
        if len(closest_points) == 0:
            self.closest_points = self.find_closest_points_to_point(
                dist_matrix)
        else:
            self.closest_points = closest_points
        self.size = init_size
        if min_r == -1:
            self.min_r = self.calc_min_r(init_size)
        else:
            min_r = min_r
        self.pvalue = pvalue
        self.dimension = dimension
        self.r = r
        self.lambda_val = lambda_val
        self.volume = volume
        self.density = density

    def remove_point(self, ind):
        if ind == self.center_point_ind:
            if len(self.closest_points) != 0:
                if self.has_center:
                    for i in range(len(self.closest_points)):
                        self.closest_points[i][1] = -1
                self.center_point_ind = self.closest_points[0][0]
                self.closest_points = self.closest_points[1:]
            self.has_center = False
            self.size -= 1
        else:
            for point in self.closest_points:
                if ind == point[0]:
                    self.closest_points.remove(point)
                    break
            self.size -= 1
        # self.volume = -1
        # self.density = -1
        # self.r = [-1, -1.0]

    def find_new_center(self, dist_matrix):
        point_ind_list = [point_ind_dist[0] for point_ind_dist in self.closest_points]
        point_ind_list.append(self.center_point_ind)
        point_ind_dist_max_dict = {point_ind: [] for point_ind in point_ind_list}
        for point_ind in point_ind_list:
            for point_ind_2 in point_ind_list:
                if point_ind == point_ind_2:
                    continue
                point_ind_dist_max_dict[point_ind].append(
                    dist_matrix.dist_matrix[point_ind, point_ind_2])
            point_ind_dist_max_dict[point_ind] = max(
                point_ind_dist_max_dict[point_ind])

        point_ind_min_max_dist = min(point_ind_dist_max_dict.keys(), key=(
            lambda k: point_ind_dist_max_dict[k]))
        if point_ind_min_max_dist != self.center_point_ind:
            self.closest_points = [point_ind_dist for point_ind_dist in self.closest_points
                                   if point_ind_dist[0] != point_ind_min_max_dist] + [[self.center_point_ind, 0.0]]
        self.center_point_ind = point_ind_min_max_dist
        self.has_center = True
        self.closest_points = self.find_closest_points_to_point(
            dist_matrix, inside_neighb=True)
        self.r = self.closest_points[-1]
        self.size = len(self.closest_points) + 1
        self.volume = self.r[1]**self.dimension
        self.density = self.size/self.volume

    def calc_min_r(self, init_size):
        min_r = 0.0
        # for i in range(init_size):
        #     min_r += self.closest_points[i][1]
        # min_r = min_r/init_size
        i = 0
        while min_r == 1E-299 or min_r == 0.0:
            min_r = self.closest_points[init_size - 1 + i][1]
            i += 1
        return min_r
    # Step 1

    def find_closest_points_to_point(self, dist_matrix, point_ind=None, inside_neighb=False):
        if point_ind == None:
            point_ind = self.center_point_ind
        i_dist_list = []
        if inside_neighb:
            #closest_points_idx_list = [i_dist[0] for i_dist in self.closest_points + [[self.center_point_ind, 0.0]]]
            i_dist_list = [[cl_point_ind[0], dist_matrix.dist_matrix[point_ind, cl_point_ind[0]]]
                           for cl_point_ind in self.closest_points + [[self.center_point_ind, 0.0]]]
        else:
            i_dist_list = [[i, dist] for i, dist in enumerate(
                dist_matrix.dist_matrix[point_ind])]
        sorted_dist = sorted(i_dist_list, key=lambda x: (
            x[1], x[0]))  # Sorting by distance
        for ind_dist in sorted_dist:
            if point_ind == ind_dist[0]:
                sorted_dist.remove(ind_dist)
                break
        return sorted_dist
    # Step 2 3

    def calc_pval_for_cur_size(self, size):
        i = size
        r_i = self.closest_points[i][1]
        r_i_next = self.closest_points[i+1][1]
        if r_i_next - self.min_r == 0.0 or r_i_next == 1E-299:
            self.pvalue_list.append(0.99)
            self.dim_i_next_list.append(0.0)
        else:
            # Hausdorff fractal dimension
            # but for two dim case it is 2
            #dim_i_next =  2
            dim_i_next = math.log(float(i) + 2.0) / \
                (math.log(2.0*r_i_next) - math.log(self.min_r))
            # Testing Poisson enrichment
            lambda_i_next = (i + 2)/(r_i_next**dim_i_next)  # Сircle area
            #pvalue = poisson.sf(i + 1, (float(i)+2.0) * (r_i/r_i_next)**dim_i_next)
            pvalue = poisson.sf(i + 1, lambda_i_next*(r_i**dim_i_next))
            if pvalue == 0.0:
                pvalue = 1E-299

            self.pvalue_list.append(pvalue)
            self.dim_i_next_list.append(dim_i_next)

    # Step 4
    # Calc min pval in 1 step
    def find_enrichment_area_radius(self, zone_size_percent=50, q_val=3):
        self.pvalue_list = []
        self.dim_i_next_list = []
        init_size = self.size
        # Search minimum pvalue at zone which set by user
        for i in range(int(len(self.closest_points)*(zone_size_percent/100)) - init_size - 1):
            self.calc_pval_for_cur_size(self.size)
            self.size += 1
            # if i > 1:
            #     mLnPvalKsd = np.std(np.array(list(map(lambda  x: -math.log10(x), pvalue_list[:-1]))))
            #     if (math.log10(pvalue_list[-2]) - math.log10(pvalue_list[-1])) > q_val*mLnPvalKsd:
            #         self.pvalue = min(pvalue_list[:-1])
            #         return self.closest_points[self.size]
        self.size += 1

        # Calc min pval in 1 step
        # Search minimum pvalue at zone which set by user
        #self.pvalue = min(self.pvalue_list[:int(len(self.pvalue_list)*(zone_size_percent/100))])
        self.pvalue = min(self.pvalue_list)
        min_pvalue_index = self.pvalue_list.index(self.pvalue)
        self.dimension = self.dim_i_next_list[min_pvalue_index]
        self.size -= len(self.pvalue_list) - min_pvalue_index - 1
        self.closest_points = self.closest_points[:self.size - 1]
        # Index of point and distance from center to it
        self.r = self.closest_points[-1]
        # Find lambda
        self.lambda_val = self.lambda_calc()
        # Find volume
        self.volume = self.r[1]**self.dimension
        self.density = self.size/self.volume

    # Calc min pval in 2 steps
    def find_enrichment_area_radius_2_steps(self, zone_size_percent=50, q_val=3, pval_thresh=0.000001):
        self.pvalue_list = []
        self.dim_i_next_list = []
        init_size = self.size
        # Search minimum pvalue at zone which set by user
        for i in range(int(len(self.closest_points)*(zone_size_percent/100)) - init_size - 1):
            self.calc_pval_for_cur_size(self.size)
            self.size += 1
            # if i > 1:
            #     mLnPvalKsd = np.std(np.array(list(map(lambda  x: -math.log10(x), pvalue_list[:-1]))))
            #     if (math.log10(pvalue_list[-2]) - math.log10(pvalue_list[-1])) > q_val*mLnPvalKsd:
            #         self.pvalue = min(pvalue_list[:-1])
            #         return self.closest_points[self.size]
        self.size += 1

        # Calc min pval in 1 step
        self.pvalue = min(self.pvalue_list)
        min_pvalue_index = self.pvalue_list.index(self.pvalue)

        # 2 step
        for i, pval in enumerate(self.pvalue_list):
            pval_thresh = 1/(i+4)
            if pval <= pval_thresh:
                self.pvalue = pval
                min_pvalue_index = i
                break
        #TODO: Debug
        # with open("debug_pval.txt", 'w') as debug_out:
        #     for i, pval in enumerate(self.pvalue_list):
        #         pval_thresh = 1/(i+4)
        #         debug_out.write(str(i) + '\t' + str(pval) + '\t' + str(pval_thresh) + '\n')
        # exit()

        # tmp_pvalue = min(self.pvalue_list[:min_pvalue_index])
        # if tmp_pvalue <= pval_thresh:
        #     self.pvalue = tmp_pvalue
        #     min_pvalue_index = self.pvalue_list.index(self.pvalue)
        import numpy
        self.dimension = self.dim_i_next_list[min_pvalue_index]
        self.size -= len(self.pvalue_list) - min_pvalue_index - 1
        self.closest_points = self.closest_points[:self.size - 1]
        # Index of point and distance from center to it
        self.r = self.closest_points[-1]
        # Find lambda
        self.lambda_val = self.lambda_calc()
        # Find volume
        self.volume = self.r[1]**self.dimension
        # Find density
        self.density = self.size/self.volume

    # Calc min pval in 2 steps
    # def find_enrichment_area_radius_2_steps(self, q_val=3):
    #     self.pvalue_list = []
    #     self.dim_i_next_list = []
    #     init_size = self.size
    #     for i in range(len(self.closest_points) - init_size - 1):
    #         self.calc_pval_for_cur_size(self.size)
    #         self.size += 1
    #         # if i > 1:
    #         #     mLnPvalKsd = np.std(np.array(list(map(lambda  x: -math.log10(x), pvalue_list[:-1]))))
    #         #     if (math.log10(pvalue_list[-2]) - math.log10(pvalue_list[-1])) > q_val*mLnPvalKsd:
    #         #         self.pvalue = min(pvalue_list[:-1])
    #         #         return self.closest_points[self.size]
    #     self.size += 1

    #     #Calc min pval in 2 steps
    #     pval_1 = min(self.pvalue_list)
    #     min_pval_1_index = self.pvalue_list.index(pval_1)
    #     if min_pval_1_index + 1 < len(self.pvalue_list):
    #         pval_2 = min(self.pvalue_list[min_pval_1_index+1:])
    #         min_pval_2_index = self.pvalue_list.index(pval_2)
    #         min_pvalue_index = -1
    #         if pval_1 == 0.0 or pval_2/pval_1 > 2:
    #             self.pvalue = pval_1
    #             min_pvalue_index = min_pval_1_index
    #         else:
    #             self.pvalue = pval_2
    #             min_pvalue_index = min_pval_2_index
    #     else:
    #         self.pvalue = pval_1
    #         min_pvalue_index = min_pval_1_index
    #     self.dimension = self.dim_i_next_list[min_pvalue_index]
    #     self.size -= len(self.pvalue_list) - min_pvalue_index - 1
    #     self.closest_points = self.closest_points[:self.size - 1]
    #     self.r = self.closest_points[-1] # Index of point and distance from center to it
    def find_pvalue_list_for_known_area(self, init_size):
        self.pvalue_list = []
        self.dim_i_next_list = []
        for i in range(len(self.closest_points) - init_size - 1):
            self.calc_pval_for_cur_size(init_size + i)

    def lambda_calc(self, lambda_init_val=0.00000000000001, r=None, size=None, pvalue=None, dim=None):
        if r == None:
            r = self.r[1]
            size = self.size
            pvalue = self.pvalue
            dim = self.dimension
        return scipy.optimize.fsolve(lambda_func, lambda_init_val, args=(r, size, pvalue, dim))[0]

    def calc_pval_for_cur_size_small_circle(self, small_circle_points_list, size):
        pvalue = self.pvalue * poisson.sf(size,
                                          self.lambda_val*(small_circle_points_list[-1][1]**self.dimension))
        if pvalue == 0.0:
            pvalue = 1E-299
        return pvalue

    def del_points_with_bad_pvalue(self, dist_matrix, small_circle_size=7):
        if self.size <= small_circle_size + 1:
            return
        #lambda_points_dict = {}
        pvalue_points_dict = {
            point_i_dist[0]: 0.99 for point_i_dist in self.closest_points}
        for point_i_dist in self.closest_points:
            small_circle_points_list = \
                self.find_closest_points_to_point(
                    dist_matrix, point_i_dist[0], inside_neighb=True)[:small_circle_size - 1]
            # Find P-value
            #min_r = calc_min_r_in_small_circle(small_circle_points_list)
            pvalue = self.calc_pval_for_cur_size_small_circle(
                small_circle_points_list, small_circle_size)
            if pvalue < pvalue_points_dict[point_i_dist[0]]:
                pvalue_points_dict[point_i_dist[0]] = pvalue
            # lambda_points_dict[point_i_dist[0]] = self.lambda_calc(lambda_init_val=0.00000000000001, \
            #                                             r=small_circle_points_list[-2][1], \
            #                                             size=small_circle_size, \
            #                                             pvalue=pvalue, \
            #                                             dim=self.dimension)
        # for point_ind in lambda_points_dict:
        #     if lambda_points_dict[point_ind] < self.lambda_val:
        #         self.remove_point(point_ind)
        for point_ind in pvalue_points_dict:
            if pvalue_points_dict[point_ind] > 0.000001:
                self.remove_point(point_ind)

    def add_points_no_filter_to_neighb(self, points_tuple, points_tuple_no_filter, points_ind_set_skip=None):
        for i, point in enumerate(points_tuple_no_filter):
            if points_ind_set_skip != None and i + len(points_tuple) in points_ind_set_skip:
                continue
            dist = points_tuple[self.center_point_ind].calc_dist(point)
            if dist < self.r[1]:
                self.closest_points.append([i + len(points_tuple), dist])
                self.size += 1
        self.closest_points = sorted(
            self.closest_points, key=lambda x: x[1])  # Sorting by distance
        # Density calculation for multi dim case
        self.density = self.size/(self.volume)

def calc_min_r_in_small_circle(small_circle_points_list):
    min_r = 1E-299
    i = 0
    while min_r == 1E-299 or min_r == 0.0:
        min_r = small_circle_points_list[i][1]
        i += 1
        if i == len(small_circle_points_list):
            return min_r
    return min_r

# def calc_pval_for_cur_size_small_circle(small_circle_points_list, size, min_r):
#     i = size - 2
#     r_i = small_circle_points_list[i][1]
#     r_i_next = small_circle_points_list[i+1][1]
#     if r_i_next - min_r == 0.0 or r_i_next == 2E-19:
#         return 0.99
#     else:
#         #Hausdorff fractal dimension
#         # but for two dim case it is 2
#         #dim_i_next =  2
#         dim_i_next = math.log(float(i) + 2.0)/(math.log(2*r_i_next) - math.log(min_r))
#         #Testing Poisson enrichment
#         lambda_i_next = (i + 2)/(r_i_next**dim_i_next) #Сircle area
#         #pvalue = poisson.sf(i + 1, (float(i)+2.0) * (r_i/r_i_next)**dim_i_next)
#         pvalue = poisson.sf(i + 1, lambda_i_next*(r_i**dim_i_next))
#         return pvalue


def lambda_func(x, *args):
    r = args[0]
    size = args[1]
    pvalue = args[2]
    dim = args[3]
    return ((size * (math.log(x) + dim*math.log(r)) -
             log_fact(size) - math.log(pvalue))/math.pow(r, dim)) - x

def log_fact(x):
    sum = 0
    for i in range(x):
        sum = sum + math.log(i+1)
    return sum

def pval_calc(lambda_val, volume, size):
    ln_pval = size*math.log(lambda_val*volume) \
        + (-lambda_val*volume) \
        - log_fact(size)
    return math.pow(math.e, ln_pval)

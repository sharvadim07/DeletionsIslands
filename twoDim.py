import math
import space
import neighborhood

class TwoDimPoint(space.Point):
    def __init__(self, x1, y1):
        self.x1 = float(x1)
        self.y1 = float(y1)

    def get_value(self):
        return (self.x1, self.y1)

    def calc_dist(self, point):
        dist = math.sqrt(((point.x1 - self.x1)**2) + ((point.y1 - self.y1)**2))
        if dist == 0.0:
            dist = 1E-299
        return dist

class TwoDimPointDel(TwoDimPoint):
    def __init__(self, x1, y1, psi1, sd1, passage_num, replicate_num, coverage, junctions_per_pos):
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.psi1 = float(psi1)
        self.sd1 = float(sd1)
        self.singPsi1 = self.psi1/self.sd1
        self.passage_num = int(passage_num)
        self.replicate_num = int(replicate_num)
        self.coverage = int(coverage)
        self.junctions_per_pos = int(junctions_per_pos)

class NeighborhoodTwoDim(neighborhood.Neighborhood):
    def init_neighborhood_from_info(self, points_tuple, radius, pvalue):
        self.r = [-1, radius]
        self.size = 0
        self.pvalue = pvalue
        self.closest_points = []
        for i, point in enumerate(points_tuple):
            if ((point.x1 - points_tuple[self.center_point_ind].x1)**2 + (point.y1 - points_tuple[self.center_point_ind].y1)**2) <= self.r[1]**2:
                self.closest_points.append(
                    [i, point.calc_dist(points_tuple[self.center_point_ind])])
                self.size += 1
        self.closest_points = sorted(self.closest_points, key=lambda x: x[1])[
            1:]  # Sorting by distance

    def calc_significance(self, points_tuple):
        # average_sign_psi = sum([points_tuple[point_ind_dist[0]].singPsi1 \
        #     for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])/float(self.size)
        average_sign_psi = 0
        sum_of_diff_sign_psi_and_avg = sum([points_tuple[point_ind_dist[0]].singPsi1 - average_sign_psi
                                            for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])
        self.sign_neighborhood = sum_of_diff_sign_psi_and_avg / \
            math.sqrt(self.size)

    def calc_sum_cov_junc_pos_per_passage_replicate(self, points_tuple, passage_num_list, replicate_num_list):
        cov_per_pass_repl_dict = {}
        junc_per_pos_per_pass_repl_dict = {}
        for passage_num in passage_num_list:
            cov_per_pass_repl_dict[passage_num] = {}
            junc_per_pos_per_pass_repl_dict[passage_num] = {}
            for replicate_num in replicate_num_list:
                sum_of_coverage = sum([points_tuple[point_ind_dist[0]].coverage
                                       if points_tuple[point_ind_dist[0]].passage_num == passage_num
                                       and points_tuple[point_ind_dist[0]].replicate_num == replicate_num else 0
                                       for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])
                cov_per_pass_repl_dict[passage_num][replicate_num] = sum_of_coverage
                sum_of_junctions_per_pos = sum([points_tuple[point_ind_dist[0]].junctions_per_pos
                                                if points_tuple[point_ind_dist[0]].passage_num == passage_num
                                                and points_tuple[point_ind_dist[0]].replicate_num == replicate_num else 0
                                                for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])
                junc_per_pos_per_pass_repl_dict[passage_num][replicate_num] = sum_of_junctions_per_pos
        return cov_per_pass_repl_dict, junc_per_pos_per_pass_repl_dict

    def calc_significance_per_passage(self,
                                      points_tuple,
                                      passage_num_list,
                                      replicate_num_list,
                                      passage_replicate_dict):
        # average_sign_psi = sum([points_tuple[point_ind_dist[0]].singPsi1 \
        #     for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])/float(self.size)
        sign_per_pass_dict = {}
        sum_sign_per_pass_dict = {}
        num_points_per_pass_dict = {}
        for passage_num in passage_num_list:
            points_num_in_cur_passage = sum([1 if points_tuple[point_ind_dist[0]].passage_num == passage_num else 0
                                             for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])
            num_points_per_pass_dict[passage_num] = points_num_in_cur_passage
            if points_num_in_cur_passage == 0:
                sign_per_pass_dict[passage_num] = 0.0
                sum_sign_per_pass_dict[passage_num] = 0.0
            else:
                average_sign_psi = 0.0
                # average_sign_psi = sum([points_tuple[point_ind_dist[0]].singPsi1 \
                #     if points_tuple[point_ind_dist[0]].passage_num == passage_num
                #     and points_tuple[point_ind_dist[0]].replicate_num == replicate_num else 0.0 \
                #     for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])/float(points_num_in_cur_passage_replicate)
                sum_of_diff_sign_psi_and_avg = sum([points_tuple[point_ind_dist[0]].singPsi1 - average_sign_psi
                                                    if points_tuple[point_ind_dist[0]].passage_num == passage_num else 0.0
                                                    for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])
                # Normalize by replicates num
                sum_of_diff_sign_psi_and_avg = sum_of_diff_sign_psi_and_avg / \
                    len(passage_replicate_dict[passage_num])
                sum_sign_per_pass_dict[passage_num] = sum_of_diff_sign_psi_and_avg
                sign_per_pass_dict[passage_num] = sum_of_diff_sign_psi_and_avg / \
                    math.sqrt(points_num_in_cur_passage)
        return sign_per_pass_dict, sum_sign_per_pass_dict, num_points_per_pass_dict

    def calc_significance_per_passage_replicate(self, points_tuple, passage_num_list, replicate_num_list):
        # average_sign_psi = sum([points_tuple[point_ind_dist[0]].singPsi1 \
        #     for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])/float(self.size)
        sign_per_pass_repl_dict = {}
        sum_sign_per_pass_repl_dict = {}
        num_points_per_pass_repl_dict = {}
        for passage_num in passage_num_list:
            sign_per_pass_repl_dict[passage_num] = {}
            sum_sign_per_pass_repl_dict[passage_num] = {}
            num_points_per_pass_repl_dict[passage_num] = {}
            for replicate_num in replicate_num_list:
                points_num_in_cur_passage_replicate = sum([1 if points_tuple[point_ind_dist[0]].passage_num == passage_num
                                                           and points_tuple[point_ind_dist[0]].replicate_num == replicate_num else 0
                                                           for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])
                num_points_per_pass_repl_dict[passage_num][replicate_num] = points_num_in_cur_passage_replicate
                if points_num_in_cur_passage_replicate == 0:
                    sign_per_pass_repl_dict[passage_num][replicate_num] = 0.0
                    sum_sign_per_pass_repl_dict[passage_num][replicate_num] = 0.0
                else:
                    average_sign_psi = 0.0
                    # average_sign_psi = sum([points_tuple[point_ind_dist[0]].singPsi1 \
                    #     if points_tuple[point_ind_dist[0]].passage_num == passage_num
                    #     and points_tuple[point_ind_dist[0]].replicate_num == replicate_num else 0.0 \
                    #     for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])/float(points_num_in_cur_passage_replicate)
                    sum_of_diff_sign_psi_and_avg = sum([points_tuple[point_ind_dist[0]].singPsi1 - average_sign_psi
                                                        if points_tuple[point_ind_dist[0]].passage_num == passage_num
                                                        and points_tuple[point_ind_dist[0]].replicate_num == replicate_num else 0.0
                                                        for point_ind_dist in self.closest_points + [[self.center_point_ind, 0.0]]])
                    sum_sign_per_pass_repl_dict[passage_num][replicate_num] = sum_of_diff_sign_psi_and_avg
                    sign_per_pass_repl_dict[passage_num][replicate_num] = sum_of_diff_sign_psi_and_avg/math.sqrt(
                        points_num_in_cur_passage_replicate)
        return sign_per_pass_repl_dict, sum_sign_per_pass_repl_dict, num_points_per_pass_repl_dict

def circle_intersect(x1, y1, x2, y2, r1, r2):
    distSq = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)
    radSumSq = (r1 + r2) * (r1 + r2)
    if (distSq == radSumSq):
        return 1
    elif (distSq > radSumSq):
        return -1
    else:
        return 0

def get_neighborhoods_without_intersection(space, Q):
    # Find neighborhood which not have intersection
    neighborhood_to_del_list = []
    for i, neighb in enumerate(space.neighb_sort_by_pval_list):
        if neighb in neighborhood_to_del_list:
            continue
        for neighb_next in space.neighb_sort_by_pval_list[i+1:]:
            if neighb_next in neighborhood_to_del_list:
                continue
            # For two dim case only
            if circle_intersect(space.points_tuple[neighb.center_point_ind].x1,
                                space.points_tuple[neighb.center_point_ind].y1,
                                space.points_tuple[neighb_next.center_point_ind].x1,
                                space.points_tuple[neighb_next.center_point_ind].y1,
                                neighb.r[1],
                                neighb_next.r[1]) == 0:
                count_intersect_points = 0
                for point_neighb_next in neighb_next.closest_points:
                    if ((space.points_tuple[point_neighb_next[0]].x1 - space.points_tuple[neighb.center_point_ind].x1)**2 +
                            (space.points_tuple[point_neighb_next[0]].y1 - space.points_tuple[neighb.center_point_ind].y1)**2) < neighb.r[1]**2:
                        count_intersect_points += 1
                if ((space.points_tuple[neighb_next.center_point_ind].x1 - space.points_tuple[neighb.center_point_ind].x1)**2 +
                        (space.points_tuple[neighb_next.center_point_ind].y1 - space.points_tuple[neighb.center_point_ind].y1)**2) < neighb.r[1]**2:
                    count_intersect_points += 1
                if count_intersect_points*100/neighb_next.size >= 20 and neighb_next not in neighborhood_to_del_list:
                    neighborhood_to_del_list.append(neighb_next)
            ####
    for neighborhood_to_del in neighborhood_to_del_list:
        space.neighb_sort_by_pval_list.remove(neighborhood_to_del)


def add_all_points_without_space(points_tuple, points_tuple_no_filter, neighb_sort_by_feature_list):
    # Add other points without filtering by frequency and length
    for neighb in neighb_sort_by_feature_list:
        neighb.add_points_no_filter_to_neighb(
            points_tuple, points_tuple_no_filter)
    return tuple(points_tuple + points_tuple_no_filter)


def find_center_ind_by_coords(points_tuple, x1, y1):
    for i, point in enumerate(points_tuple):
        if point.x1 == x1 and point.y1 == y1:
            return i


def read_neighb_info_file(neighb_info_file_name):
    neighb_info_list = []
    with open(neighb_info_file_name, 'r') as neighb_info_file:
        header = []
        for i, line in enumerate(neighb_info_file):
            line_list = line.split('\t')
            if i != 0:
                if line_list[header.index('NCenterX')] != ' ':
                    neighb_info_list.append(line_list[:5])
            else:
                header = line_list
    return neighb_info_list


def init_neighb_list(neighb_info_file_name, points_tuple, init_size=3):
    neighb_info_list = read_neighb_info_file(neighb_info_file_name)
    neighb_list = []
    for neighb_info in neighb_info_list:
        x1 = float(neighb_info[0])
        y1 = float(neighb_info[1])
        radius = float(neighb_info[2])
        pvalue = float(neighb_info[4])
        neighb = NeighborhoodTwoDim(
            find_center_ind_by_coords(points_tuple, x1, y1))
        neighb.init_neighborhood_from_info(points_tuple, radius, pvalue)
        neighb_list.append(neighb)
    for neighb in neighb_list:
        neighb.find_pvalue_list_for_known_area(init_size)
    return neighb_list

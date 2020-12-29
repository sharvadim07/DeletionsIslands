import space
import neighborhood
import math
from scipy.spatial import distance

class MultiDimPoint(space.Point):
    def __init__(self, init_param_tuple, special_features_idx_dict, ln_param=False):
        self.special_features_values_dict = {}
        special_features_idx_list = [special_features_idx_dict[special_feature]
                                     for special_feature in special_features_idx_dict]
        if ln_param:
            self.param_tuple = tuple([math.log(1 + float(init_param)) for i, init_param in enumerate(init_param_tuple)
                                      if i not in special_features_idx_list])
        else:
            self.param_tuple = tuple([float(init_param) for i, init_param in enumerate(init_param_tuple)
                                      if i not in special_features_idx_list])
        for special_feature in special_features_idx_dict:
            self.special_features_values_dict[special_feature] = init_param_tuple[special_features_idx_dict[special_feature]]

    def get_value(self):
        return (self.param_tuple)

    def calc_dist(self, point):
        dist = distance.euclidean(self.param_tuple, point.param_tuple)
        if dist == 0.0:
            dist = 1E-299
        return dist

class NeighborhoodMultiDim(neighborhood.Neighborhood):
    def calc_neighborhood_avg_spec_feature(self, points_tuple, spec_feature):
        avg = float(
            points_tuple[self.center_point_ind].special_features_values_dict[spec_feature])
        for point in self.closest_points + [[self.center_point_ind, 0]]:
            avg += float(points_tuple[point[0]
                                      ].special_features_values_dict[spec_feature])
        return avg/self.size

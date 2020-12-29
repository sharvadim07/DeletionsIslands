import space
import multiDim

class ParTable:
    def __init__(self, par_table_file_name, special_features_cols, add_features_name_val_list):
        self.header = None
        self.param_tuple_list = []
        self.special_features_cols_idx_dict = {}
        with open(par_table_file_name, 'r') as par_table_file:
            for i, line in enumerate(par_table_file):
                if i == 0:
                    self.header = tuple(line.strip().split('\t') + [add_feature_name_val[0]
                                                                    for add_feature_name_val in add_features_name_val_list])
                    for special_feature in special_features_cols:
                        self.special_features_cols_idx_dict[special_feature] = \
                            self.header.index(special_feature)
                else:
                    splitted_line_tuple = tuple(line.strip().split('\t') + [add_feature_name_val[1]
                                                                            for add_feature_name_val in add_features_name_val_list])
                    if len(splitted_line_tuple) != len(self.header):
                        raise IOError(
                            'Num of elements in line not corresponds with header size!')
                    self.param_tuple_list.append(splitted_line_tuple)

    def append(self, par_table_file_name, add_features_name_val_list):
        with open(par_table_file_name, 'r') as par_table_file:
            for i, line in enumerate(par_table_file):
                if i != 0:
                    splitted_line_tuple = tuple(line.strip().split('\t') + [add_feature_name_val[1]
                                                                            for add_feature_name_val in add_features_name_val_list])
                    if len(splitted_line_tuple) != len(self.header):
                        raise IOError(
                            'Num of elements in line not corresponds with header size!')
                    self.param_tuple_list.append(splitted_line_tuple)

    def __len__(self):
        return len(self.param_tuple_list)

    def __getitem__(self, position):
        return self.param_tuple_list[position]

def generate_points_tuple(par_table):
    points_tuple = tuple([multiDim.MultiDimPoint(param_tuple, par_table.special_features_cols_idx_dict)
                          for param_tuple in par_table.param_tuple_list])
    return points_tuple

# Last run params "-K3", "-A500", "-Z30"

def param_table_processing(args, out_dir, type_of_table):
    import re
    import os

    if type_of_table == 'phydyn':
        spec_features = ['mcmc_step', 'glob_iter',
                         'residual', 'treeLikelihood']
    elif type_of_table == 'expdyn':
        spec_features = ['mcmc_step', 'glob_iter', 'residual']
    elif type_of_table == 'single_cell':
        spec_features = ['DataType', 'Cluster', 'BarCode', 'Cell_Id']

    if args.par_tables != None:
        par_table = None
        for i, par_table_file_name in enumerate(args.par_tables.split(',')):
            if par_table_file_name == '':
                break
            if type_of_table == 'phydyn' or type_of_table == 'expdyn':
                global_iteration_num = int(
                    re.search(r'_([0-9]+)', os.path.basename(par_table_file_name)).groups()[0])
                add_features_name_val_list = [
                    ['glob_iter', str(global_iteration_num)]]
            elif 'single_cell':
                #sample_id = int(re.search(r'_([0-9]+)', os.path.basename(par_table_file_name)).groups()[0])
                #add_features_name_val_list = [['sample_id', str(sample_id)]]
                add_features_name_val_list = []
            if i == 0:
                par_table = ParTable(par_table_file_name,
                                     spec_features, add_features_name_val_list)
            else:
                par_table.append(par_table_file_name,
                                 add_features_name_val_list)
        # Determine centers of neighborhoods
        points_tuple = generate_points_tuple(par_table)
        filter_points_tuple = tuple(list(points_tuple)[::10])
        #new_space = space.Space(points_tuple, args.num_K_points, args.A_value, args.zone_size_percent)
        if args.do_two_step:
            new_space = space.Space(points_tuple=filter_points_tuple,
                                    num_K_points=args.num_K_points,
                                    A=args.A_value,
                                    sort_neighb_by_density=True,
                                    neighb_post_process=False,
                                    zone_size_percent=100,
                                    num_of_random_centers_from_isalnd=1,
                                    second_step=True,
                                    min_neighb_part_from_all=0.0025)
            # Find best pvalues for point located in found neighborhoods
            new_space.find_pvalue_for_points_in_all_neighb()
            import numpy
            pvalue_points_in_neighb_med = \
                numpy.median(
                    list(new_space.pvalue_points_in_neighb_dict.values()))
            # Filter points by pvalue
            new_space.filter_points_by_pvalue(
                pvalue_threshold=pvalue_points_in_neighb_med)
            os.chdir(out_dir)
            # Debug
            # print_neighborhoods_all_points('param_neighborhoods_1_single_cell.txt',
            #                                 new_space,
            #                                 [],
            #                                 spec_features,
            #                                 par_table.header)

            points_tuple_in_cur_neighbs = space.get_points_tuple_in_cur_neghbs(
                new_space)
            points_tuple_not_in_cur_neighbs = tuple([point for point in points_tuple
                                                     if point not in points_tuple_in_cur_neighbs])

            two_step_space = space.Space(points_tuple=points_tuple_in_cur_neighbs,
                                         num_K_points=args.num_K_points,
                                         A=args.A_value,
                                         zone_size_percent=100,
                                         neighb_post_process=False,
                                         second_step=True,
                                         sort_neighb_by_density=True,
                                         num_of_random_centers_from_isalnd=1,
                                         min_neighb_part_from_all=0.0025)
            two_step_space.find_new_centers_of_neighborhoods()

            space.add_all_points(two_step_space.points_tuple,
                                 points_tuple_not_in_cur_neighbs, two_step_space)

            two_step_space.get_neighborhoods_without_intersection_by_density()

            new_space = two_step_space
        elif args.do_mod_step:
            new_space = space.Space(points_tuple=filter_points_tuple,
                                    num_K_points=args.num_K_points,
                                    A=args.A_value,
                                    sort_neighb_by_density=True,
                                    neighb_post_process=False,
                                    zone_size_percent=100,
                                    num_of_random_centers_from_isalnd=1,
                                    second_step=True,
                                    min_neighb_part_from_all=0.0025)
            # Find down quartile of lambda values of neighborhoods
            lambda_val_list = sorted(
                [neighb.lambda_val for neighb in new_space.neighb_sort_by_feature_list])
            quartile_idx = int((len(lambda_val_list)-1)/4)
            down_quartile_lambda = lambda_val_list[quartile_idx]
            
            # Reduce number of redundant neighborhoods
            new_space.get_neighborhoods_without_intersection_by_density(only_full_intersection=True)
            
            # Find best pvalues for point located in found neighborhoods
            new_space.find_pvalue_for_points_in_all_neighb()
            #import numpy
            # pvalue_points_in_neighb_med = \
            #    numpy.median(list(new_space.pvalue_points_in_neighb_dict.values()))
            p_val_list = list(new_space.pvalue_points_in_neighb_dict.values())
            quartile_idx = int((len(p_val_list)-1)/4)
            pvalue_points_up_quartile = sorted(p_val_list, reverse=True)[
                quartile_idx]
            # Filter points by pvalue
            new_space.filter_points_by_pvalue(
                pvalue_threshold=pvalue_points_up_quartile)

            points_tuple_in_cur_neighbs = space.get_points_tuple_in_cur_neghbs(
                new_space)
            points_tuple_not_in_cur_neighbs = tuple([point for point in points_tuple
                                                     if point not in points_tuple_in_cur_neighbs])

            # Find graphs of small circles
            new_space.find_small_circles_graphs(points_tuple_not_in_cur_neighbs,
                                                down_quartile_lambda)
        else:
            new_space = space.Space(points_tuple=points_tuple,
                                    num_K_points=args.num_K_points,
                                    A=args.A_value,
                                    sort_neighb_by_density=True,
                                    neighb_post_process=True,
                                    zone_size_percent=100,
                                    num_of_random_centers_from_isalnd=1,
                                    second_step=True,
                                    min_neighb_part_from_all=0.0025)

        os.chdir(out_dir)
        if type_of_table == 'phydyn':
            print_neighborhoods_all_points('param_dyn_neighborhoods_phy.txt',
                                           new_space,
                                           ['residual', 'treeLikelihood'],
                                           spec_features,
                                           par_table.header)
        elif type_of_table == 'expdyn':
            print_neighborhoods_all_points('param_dyn_neighborhoods_exp.txt',
                                           new_space,
                                           ['residual'],
                                           spec_features,
                                           par_table.header)
        elif type_of_table == 'single_cell':
            print_neighborhoods_all_points('param_neighborhoods_single_cell.txt',
                                           new_space,
                                           [],
                                           spec_features,
                                           par_table.header)

###
# Printing results
def print_neighborhoods_all_points(out_file_name,
                                   cur_space,
                                   avg_spec_features,
                                   spec_features,
                                   table_header):
    import os
    with open(os.path.basename(out_file_name), 'w') as out_file:
        for i, neighb in enumerate(cur_space.neighb_sort_by_feature_list):
            avg_info = ''
            avg_spec_features_val_dict = {}
            for avg_spec_feature in avg_spec_features:
                avg_spec_feature_val = neighb.calc_neighborhood_avg_spec_feature(
                    cur_space.points_tuple, avg_spec_feature)
                avg_info += ' avg_' + avg_spec_feature + \
                    '=' + str(avg_spec_feature_val)
                avg_spec_features_val_dict[avg_spec_feature] = avg_spec_feature_val

            out_file.write('#neighb ' + str(i) + ' info: pValue=' + str(neighb.pvalue) +
                           ' dimension=' + str(neighb.dimension) +
                           ' volume=' + str(neighb.volume) +
                           ' size=' + str(neighb.size) +
                           str(avg_info) + '\n')
            header = [col_head for col_head in table_header if col_head not in spec_features] + \
                spec_features + \
                ['neighb_id', 'neighb_pValue', 'neighb_dimension', 'neighb_volume', 'neighb_size'] + \
                ['avg_' + avg_spec_feature for avg_spec_feature in avg_spec_features]
            out_file.write('\t'.join(header) + '\n')
            for point_ind_dist in neighb.closest_points[:neighb.size] + [[neighb.center_point_ind, 0]]:
                point = cur_space.points_tuple[point_ind_dist[0]]
                line = list(map(str, point.param_tuple)) + \
                    [str(point.special_features_values_dict[spec_feature]) for spec_feature in spec_features] + \
                    [str(i), str(neighb.pvalue), str(neighb.dimension), str(neighb.volume), str(neighb.size)] + \
                    [str(avg_spec_features_val_dict[avg_spec_feature])
                        for avg_spec_feature in avg_spec_features]
                out_file.write('\t'.join(map(str, line)) + '\n')
    #space.print_dist_mat(cur_space.dist_matrix, 'dist_matrix.txt')
    # Transform dist matrix
    #trans_dist_matrix = space.transform_dist_matrix(cur_space.dist_matrix, cur_space.neighb_sort_by_feature_list)
    #space.print_dist_mat(trans_dist_matrix, 'transformed_dist_matrix.txt')
###

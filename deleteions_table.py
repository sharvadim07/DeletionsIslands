import twoDim
import space


class DelTable:
    def __init__(self,
                 del_table_file_name,
                 freq_threshold,
                 freq_threshold_2,
                 del_length,
                 passage_num,
                 replicate_num,
                 passage_replicate_dict):
        self.header = []
        self.deletions_list = []
        self.deletions_list_no_filter = []
        if freq_threshold == 'auto':
            self.freq_threshold = 0.0
            self.freq_threshold_2 = 0.0
        else:
            self.freq_threshold = float(freq_threshold)
            self.freq_threshold_2 = float(freq_threshold_2)
        self.del_length = int(del_length)
        self.passage_num_list = [passage_num]
        self.replicate_num_list = [replicate_num]
        self.passage_replicate_dict = passage_replicate_dict
        with open(del_table_file_name, 'r') as del_table_file:
            psi_col_num = -1
            pos_start_col_num = -1
            pos_end_col_num = -1
            for i, line in enumerate(del_table_file):
                if i == 0:
                    self.header = line.strip().split('\t')
                    psi_col_num = self.header.index('psi1')
                    pos_start_col_num = self.header.index('posStart')
                    pos_end_col_num = self.header.index('PosEnd')
                    self.header.append('PassageNum')
                    self.header.append('ReplicateNum')
                else:
                    line_el_list = line.strip().split('\t')
                    line_el_list.append(passage_num)
                    line_el_list.append(replicate_num)
                    if len(line_el_list) == len(self.header):
                        if float(line_el_list[psi_col_num]) >= self.freq_threshold and \
                                abs(int(line_el_list[pos_end_col_num]) - int(line_el_list[pos_start_col_num])) >= self.del_length:
                            self.deletions_list.append(line_el_list)
                        elif float(line_el_list[psi_col_num]) >= self.freq_threshold_2 and \
                                abs(int(line_el_list[pos_end_col_num]) - int(line_el_list[pos_start_col_num])) >= self.del_length:
                            self.deletions_list_no_filter.append(line_el_list)
                    else:
                        raise IOError(
                            'Num of elements in line not corresponds with header size!')

    def append_deletions_from_file(self, del_table_file_name, passage_num, replicate_num):
        with open(del_table_file_name, 'r') as del_table_file:
            psi_col_num = self.header.index('psi1')
            pos_start_col_num = self.header.index('posStart')
            pos_end_col_num = self.header.index('PosEnd')
            if passage_num not in self.passage_num_list:
                self.passage_num_list.append(passage_num)
            if replicate_num not in self.replicate_num_list:
                self.replicate_num_list.append(replicate_num)
            self.passage_num_list.sort()
            self.replicate_num_list.sort()
            for i, line in enumerate(del_table_file):
                if i != 0:
                    line_el_list = line.strip().split('\t')
                    line_el_list.append(passage_num)
                    line_el_list.append(replicate_num)
                    if len(line_el_list) == len(self.header):
                        if float(line_el_list[psi_col_num]) >= self.freq_threshold and \
                                abs(int(line_el_list[pos_end_col_num]) - int(line_el_list[pos_start_col_num])) >= self.del_length:
                            self.deletions_list.append(line_el_list)
                        elif float(line_el_list[psi_col_num]) >= self.freq_threshold_2 and \
                                abs(int(line_el_list[pos_end_col_num]) - int(line_el_list[pos_start_col_num])) >= self.del_length:
                            self.deletions_list_no_filter.append(line_el_list)
                    else:
                        raise IOError(
                            'Num of elements in line not corresponds with header size!')

    def generate_points_tuple(self, deletions_list):
        if len(self.deletions_list) > 10:
            return tuple(twoDim.TwoDimPointDel(deletion[self.header.index('posStart')],
                                               deletion[self.header.index(
                                                   'PosEnd')],
                                               deletion[self.header.index(
                                                   'psi1')],
                                               deletion[self.header.index(
                                                   'sd1')],
                                               deletion[self.header.index(
                                                   'PassageNum')],
                                               deletion[self.header.index(
                                                   'ReplicateNum')],
                                               deletion[self.header.index(
                                                   'Coverage')],
                                               deletion[self.header.index('JunctionsPerPos')]) for deletion in deletions_list
                         if float(deletion[self.header.index('Coverage')]) > 0)
        else:
            raise ValueError('deletions_list is too small for working!')


def check_point_coords_in_list(c_point, points_list):
    for point in points_list:
        if point.x1 == c_point.x1 and point.y1 == c_point.y1:
            return True
    return False


def get_not_duplicated_points_tuple(points_tuple):
    points_no_duplicated_list = []
    for point in points_tuple:
        if not check_point_coords_in_list(point, points_no_duplicated_list):
            points_no_duplicated_list.append(point)
    return tuple(points_no_duplicated_list)


def deletions_table_processing(args, out_dir):
    import re
    import os
    import shutil
    import numpy
    if args.del_tables != None and args.freq_filt != None and args.del_length != None:
        del_table_obj = None
        passage_replicate_dict = {}
        for i, del_table in enumerate(args.del_tables.split(',')):
            if del_table == '':
                break
            #passage_num = int(re.search(r'_p([0-9]+)_', os.path.basename(del_table)).groups()[0])
            #replicate_num = int(re.search(r'_H_([0-9]+)', os.path.basename(del_table)).groups()[0])
            passage_num = int(
                re.search(r'group_([0-9]+)_', os.path.basename(del_table)).groups()[0])
            replicate_num = int(
                re.search(r'sample_([0-9]+)_', os.path.basename(del_table)).groups()[0])
            if passage_num not in passage_replicate_dict:
                passage_replicate_dict[passage_num] = [replicate_num]
            else:
                passage_replicate_dict[passage_num].append(replicate_num)

            if i == 0:
                del_table_obj = DelTable(del_table,
                                         args.freq_filt,
                                         args.freq_filt_2,
                                         args.del_length,
                                         passage_num,
                                         replicate_num,
                                         passage_replicate_dict)
            else:
                del_table_obj.append_deletions_from_file(
                    del_table, passage_num, replicate_num)
        # Determine centers of neighborhoods
        points_tuple = del_table_obj.generate_points_tuple(
            del_table_obj.deletions_list)
        points_tuple_no_filter = del_table_obj.generate_points_tuple(
            del_table_obj.deletions_list_no_filter)
        freq_filt = args.freq_filt
        # Auto setting of deletions frequency
        if args.freq_filt == 'auto':
            psi1_list = [point.psi1 for point in points_tuple]
            median_frq = numpy.median(psi1_list)
            while len(psi1_list) > 10000 and median_frq < 0.9:
                median_frq = median_frq*2
                psi1_list = [psi1 for psi1 in psi1_list if psi1 >= median_frq]
            points_tuple_no_filter = tuple(
                point for point in points_tuple if point.psi1 < median_frq)
            points_tuple = tuple(
                point for point in points_tuple if point.psi1 >= median_frq)
            freq_filt = median_frq
            del_table_obj.freq_threshold = freq_filt

        # Write to out
        # shutil.rmtree(out_dir)
        # os.mkdir(out_dir)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        os.chdir(out_dir)
        if args.neighb_info != None:
            # TODO: don't tested with last changes of program
            if args.only_coords:
                neighborhood_list = twoDim.init_neighb_list(
                    args.neighb_info, points_tuple)
                print_neighborhoods_points_for_constr_coords('constructed_neighborhoods.txt',
                                                             neighborhood_list,
                                                             points_tuple)
                print_all_points(
                    'all_points_frq_' + str(del_table_obj.freq_threshold), points_tuple)
                new_points_tuple = twoDim.add_all_points_without_space(points_tuple,
                                                                       points_tuple_no_filter,
                                                                       neighborhood_list)
                print_neighborhoods_all_points('coords_points',
                                               neighborhood_list,
                                               new_points_tuple)
                print_all_points(
                    'all_points_frq_' + str(del_table_obj.freq_threshold_2), new_points_tuple)
            else:
                new_points_tuple = tuple(points_tuple + points_tuple_no_filter)
                neighborhood_list = twoDim.init_neighb_list(
                    args.neighb_info, new_points_tuple)
                print_neighborhoods_info('neighb_info.txt',
                                         neighborhood_list,
                                         new_points_tuple,
                                         del_table_obj.passage_num_list,
                                         del_table_obj.replicate_num_list,
                                         del_table_obj.passage_replicate_dict)
        else:
            print_all_points('all_points_frq_' +
                             str(del_table_obj.freq_threshold), points_tuple)
            # Second step neighb identification
            if args.do_two_step:
                # Use only not duplicated points by coords
                points_no_duplicated_tuple = get_not_duplicated_points_tuple(
                    points_tuple)
                print_all_points('all_points_not_dupl_frq_' +
                                 str(del_table_obj.freq_threshold), points_no_duplicated_tuple)
                new_space = space.Space(points_tuple=points_no_duplicated_tuple,
                                        num_K_points=args.num_K_points,
                                        A=args.A_value,
                                        sort_neighb_by_density=True,
                                        neighb_post_process=False,
                                        zone_size_percent=100,
                                        num_of_random_centers_from_isalnd=1,
                                        second_step=True,
                                        min_neighb_part_from_all=0.0025)
                print_neighborhoods_points_for_constr_coords('constructed_neighborhoods_step_1.txt',
                                                             new_space.neighb_sort_by_feature_list,
                                                             new_space.points_tuple)

                # Detection of the most enriched neighborhood for each center
                # Resolving intersection in two dim case (not used now)
                # twoDim.get_neighborhoods_without_intersection(new_space, 1000)

                # Resolving intersection in all cases (multi dim and two dim)
                # Find intersection and delete points from neighborhood with smaller denisty
                # new_space.get_neighborhoods_without_intersection_by_density()

                # Find best pvalues for point located in found neighborhoods
                new_space.find_pvalue_for_points_in_all_neighb()
                pvalue_points_in_neighb_med = \
                    numpy.median(
                        list(new_space.pvalue_points_in_neighb_dict.values()))
                # Filter points by pvalue
                new_space.filter_points_by_pvalue(
                    pvalue_threshold=pvalue_points_in_neighb_med)

                print_neighborhoods_points_for_constr_coords_with_pval('constructed_clusters_step_1_filter_by_pval.txt',
                                                                       new_space.neighb_sort_by_feature_list,
                                                                       new_space.points_tuple,
                                                                       new_space.pvalue_points_in_neighb_dict)

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
                                             min_neighb_part_from_all=0.0025)

                print_neighborhoods_points_for_constr_coords('constructed_neighborhoods_step_2.txt',
                                                             two_step_space.neighb_sort_by_feature_list,
                                                             two_step_space.points_tuple)
                two_step_space.find_new_centers_of_neighborhoods()
                print_neighborhoods_points_for_constr_coords('constructed_neighborhoods_new_centers_step_2.txt',
                                                             two_step_space.neighb_sort_by_feature_list,
                                                             two_step_space.points_tuple)
                # After neighborhood construction we use all points
                # for next steps of clustering
                new_points_tuple = space.add_all_points(points_tuple_in_cur_neighbs,
                                                        points_tuple_not_in_cur_neighbs,
                                                        two_step_space)

                if len(points_tuple_no_filter) > 0:
                    new_points_tuple = space.add_all_points(new_points_tuple,
                                                            points_tuple_no_filter,
                                                            two_step_space)
                    two_step_space.points_tuple = new_points_tuple
                    print_all_points(
                        'all_points_frq_' + str(del_table_obj.freq_threshold_2), new_points_tuple)
                    print_neighborhoods_all_points('coords_points_step_2',
                                                   two_step_space.neighb_sort_by_feature_list,
                                                   two_step_space.points_tuple)
                # Resolving intersection in all cases (multi dim and two dim)
                # Find intersection and delete points from neighborhood with smaller denisty
                two_step_space.get_neighborhoods_without_intersection_by_density()
                print_neighborhoods_points_for_constr_coords('constructed_clusters_step_2_with_all_points.txt',
                                                             two_step_space.neighb_sort_by_feature_list,
                                                             two_step_space.points_tuple)

                new_space = two_step_space
            elif args.do_mod_step:
                points_no_duplicated_tuple = get_not_duplicated_points_tuple(
                    points_tuple)
                print_all_points('all_points_not_dupl_frq_' +
                                 str(del_table_obj.freq_threshold), points_no_duplicated_tuple)
                new_space = space.Space(points_tuple=points_no_duplicated_tuple,
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
                #import numpy
                # pvalue_points_in_neighb_med = \
                #    numpy.median(list(new_space.pvalue_points_in_neighb_dict.values()))
                p_val_list = list(
                    new_space.pvalue_points_in_neighb_dict.values())
                pvalue_points_up_quartile = sorted(p_val_list, reverse=True)[
                    int(len(p_val_list)/4) + 1]
                # Filter points by pvalue
                new_space.filter_points_by_pvalue(
                    pvalue_threshold=pvalue_points_up_quartile)

                points_tuple_in_cur_neighbs = space.get_points_tuple_in_cur_neghbs(
                    new_space)
                points_tuple_not_in_cur_neighbs = tuple([point for point in tuple(points_tuple + points_tuple_no_filter)
                                                         if point not in points_tuple_in_cur_neighbs])

                # Find graphs of small circles
                new_space.find_small_circles_graphs(
                    points_tuple_not_in_cur_neighbs)
                print_neighborhoods_points_for_constr_coords('constructed_clusters_step_2_with_all_points.txt',
                                                             new_space.neighb_sort_by_feature_list,
                                                             new_space.points_tuple)
            else:
                new_space = space.Space(points_tuple=points_no_duplicated_tuple,
                                        num_K_points=args.num_K_points,
                                        A=args.A_value,
                                        sort_neighb_by_density=True,
                                        neighb_post_process=True,
                                        zone_size_percent=100,
                                        num_of_random_centers_from_isalnd=1,
                                        second_step=True,
                                        min_neighb_part_from_all=0.0025)
                if len(points_tuple_no_filter) > 0:
                    # Adding other points with lower freq threshold to neighb circle
                    # Check of belonging of point to neighb circle
                    new_points_tuple = space.add_all_points(points_tuple,
                                                            points_tuple_no_filter,
                                                            new_space)
                    print_all_points(
                        'all_points_frq_' + str(del_table_obj.freq_threshold_2), new_points_tuple)
                    print_neighborhoods_all_points('coords_points',
                                                   new_space.neighb_sort_by_feature_list,
                                                   new_points_tuple)
                else:
                    new_points_tuple = points_tuple

            print_neighborhoods_info('out_neighb_info.txt',
                                     new_space.neighb_sort_by_feature_list,
                                     new_space.points_tuple,
                                     del_table_obj.passage_num_list,
                                     del_table_obj.replicate_num_list,
                                     del_table_obj.passage_replicate_dict)
            # TODO: Working with dist matrix
            #new_space.dist_matrix = space.DistanceMatrix(new_points_tuple)
            #space.print_dist_mat(new_space.dist_matrix, 'dist_matrix.txt')
            # Transform dist matrix
            #trans_dist_matrix = space.transform_dist_matrix(new_space.dist_matrix, new_space.neighb_sort_by_feature_list)
            #space.print_dist_mat(trans_dist_matrix, 'transformed_dist_matrix.txt')

###
# Printing results


def print_neighborhoods_points_for_constr_coords_with_pval(out_file_name,
                                                           neighborhood_list,
                                                           points_tuple,
                                                           pvalue_points_in_neighb_dict,
                                                           init_size=3):
    import os
    with open(os.path.basename(out_file_name), 'w') as out_file:
        out_file.write(
            'Nid\tNCenterX\tNCenterY\tRadius\tSize\tpValue\tDensity\t\tCoordX\tCoordY\tDistToCenter\tpvalue_in_constr\tDimension\tpoint_pvalue\n')
        for i, neighb in enumerate(neighborhood_list):
            out_file.write(str(i) + '\t' +
                           str(points_tuple[neighb.center_point_ind].x1) + '\t' +
                           str(points_tuple[neighb.center_point_ind].y1) + '\t' +
                           str(neighb.r[1]) + '\t' +
                           str(neighb.size) + '\t' +
                           str(neighb.pvalue) + '\t' +
                           str(neighb.density) + '\t' +
                           '\t' +
                           str(points_tuple[neighb.center_point_ind].x1) + '\t' +
                           str(points_tuple[neighb.center_point_ind].y1) + '\t' +
                           '0.0' + '\t' +
                           '-' + '\t' +
                           '-' + '\t' +
                           str(pvalue_points_in_neighb_dict[neighb.center_point_ind]) + '\n')
            for i, point_ind_dist in enumerate(neighb.closest_points):
                if i < 3 or i == len(neighb.closest_points) - 1:
                    out_file.write(' \t \t \t \t \t \t \t\t' + str(points_tuple[point_ind_dist[0]].x1) + '\t' +
                                   str(points_tuple[point_ind_dist[0]].y1) + '\t' +
                                   str(point_ind_dist[1]) + '\t' +
                                   '-' + '\t' +
                                   '-' + '\t' +
                                   str(pvalue_points_in_neighb_dict[point_ind_dist[0]]) + '\n')
                else:
                    out_file.write(' \t \t \t \t \t \t \t\t' + str(points_tuple[point_ind_dist[0]].x1) + '\t' +
                                   str(points_tuple[point_ind_dist[0]].y1) + '\t' +
                                   str(point_ind_dist[1]) + '\t' +
                                   str(neighb.pvalue_list[i - init_size]) + '\t' +
                                   str(neighb.dim_i_next_list[i - init_size]) + '\t' +
                                   str(pvalue_points_in_neighb_dict[point_ind_dist[0]]) + '\n')


def print_neighborhoods_points_for_constr_coords(out_file_name,
                                                 neighborhood_list,
                                                 points_tuple,
                                                 init_size=3):
    import os
    with open(os.path.basename(out_file_name), 'w') as out_file:
        out_file.write(
            'Nid\tNCenterX\tNCenterY\tRadius\tSize\tpValue\tDensity\t\tCoordX\tCoordY\tDistToCenter\n')
        for i, neighb in enumerate(neighborhood_list):
            out_file.write(str(i) + '\t' +
                           str(points_tuple[neighb.center_point_ind].x1) + '\t' +
                           str(points_tuple[neighb.center_point_ind].y1) + '\t' +
                           str(neighb.r[1]) + '\t' +
                           str(neighb.size) + '\t' +
                           str(neighb.pvalue) + '\t' +
                           str(neighb.density) + '\t' +
                           '\t' +
                           str(points_tuple[neighb.center_point_ind].x1) + '\t' +
                           str(points_tuple[neighb.center_point_ind].y1) + '\t' +
                           '0.0' + '\n')
            for i, point_ind_dist in enumerate(neighb.closest_points):
                if i < 3 or i == len(neighb.closest_points) - 1:
                    out_file.write(' \t \t \t \t \t \t \t\t' + str(points_tuple[point_ind_dist[0]].x1) + '\t' +
                                   str(points_tuple[point_ind_dist[0]].y1) + '\t' +
                                   str(point_ind_dist[1]) + '\n')
                else:
                    out_file.write(' \t \t \t \t \t \t \t\t' + str(points_tuple[point_ind_dist[0]].x1) + '\t' +
                                   str(points_tuple[point_ind_dist[0]].y1) + '\t' +
                                   str(point_ind_dist[1]) + '\n')


def print_neighborhoods_all_points(out_file_name,
                                   neighborhood_list,
                                   points_tuple):
    import os
    for i, neighb in enumerate(neighborhood_list):
        with open(str(i) + '_' + str(points_tuple[neighb.center_point_ind].x1) + '_' +
                  str(points_tuple[neighb.center_point_ind].y1) + '_'
                  + os.path.basename(out_file_name) + '.txt', 'w') as out_file:
            out_file.write(str(points_tuple[neighb.center_point_ind].x1) + '\t' +
                           str(points_tuple[neighb.center_point_ind].y1) + '\n')
            for point_ind_dist in neighb.closest_points:
                out_file.write(str(points_tuple[point_ind_dist[0]].x1) + '\t' +
                               str(points_tuple[point_ind_dist[0]].y1) + '\n')


def print_all_points(out_file_name, points_tuple):
    import os
    with open(out_file_name + '.txt', 'w') as out_file:
        for point in points_tuple:
            out_file.write(str(point.x1) + '\t' + str(point.y1) + '\n')


def print_neighborhoods_info(out_file_name,
                             neighb_sort_by_feature_list,
                             points_tuple,
                             passage_num_list,
                             replicate_num_list,
                             passage_replicate_dict):
    import os
    with open('significance_' + os.path.basename(out_file_name), 'w') as out_file1, \
            open('sum_significance_' + os.path.basename(out_file_name), 'w') as out_file2, \
            open('num_points_' + os.path.basename(out_file_name), 'w') as out_file3, \
            open('cov_junct_per_pos_' + os.path.basename(out_file_name), 'w') as out_file4:
        sign_passage_header_list = [
            'signPass' + str(passage) for passage in passage_num_list]
        sign_sum_passage_header_list = [
            'signSumPass' + str(passage) for passage in passage_num_list]
        num_points_passage_header_list = [
            'numPointsPass' + str(passage) for passage in passage_num_list]
        coverage_passage_header_list = [
            'coveragePass' + str(passage) for passage in passage_num_list]
        junct_per_pos_passage_header_list = [
            'junctPerPosPass' + str(passage) for passage in passage_num_list]
        out_file1.write('NCenterX\tNCenterY\tRadius\tSize\tpValue\t' +
                        '\t'.join(map(str, sign_passage_header_list)) + '\n')
        out_file2.write('NCenterX\tNCenterY\tRadius\tSize\tpValue\t' +
                        '\t'.join(map(str, sign_sum_passage_header_list)) + '\n')
        out_file3.write('NCenterX\tNCenterY\tRadius\tSize\tpValue\t' +
                        '\t'.join(map(str, num_points_passage_header_list)) + '\n')
        out_file4.write('NCenterX\tNCenterY\tRadius\tSize\tpValue\t' +
                        '\t'.join(map(str, coverage_passage_header_list)) + '\t' +
                        '\t'.join(map(str, junct_per_pos_passage_header_list)) + '\n')
        for neighb in neighb_sort_by_feature_list:
            if neighb.r != 0.0:
                # Calculate significance
                # Per passage
                sign_per_pass_dict, \
                    sum_sign_per_pass_dict, \
                    num_points_per_pass_dict = \
                    neighb.calc_significance_per_passage(points_tuple,
                                                         passage_num_list,
                                                         replicate_num_list,
                                                         passage_replicate_dict)
                # Per passage replicate
                # sign_per_pass_repl_dict, \
                # sum_sign_per_pass_repl_dict, \
                # num_points_per_pass_repl_dict = \
                #     neighb.calc_significance_per_passage_replicate(points_tuple, passage_num_list, replicate_num_list)
                # Calculate coverage and junctions per position
                cov_per_pass_repl_dict, \
                    junc_per_pos_per_pass_repl_dict = \
                    neighb.calc_sum_cov_junc_pos_per_passage_replicate(
                        points_tuple, passage_num_list, replicate_num_list)
                # Calc values coverage and junctions per position per passage
                list_of_passage_cov = []
                list_of_passage_junc_per_pos = []
                for passage_num in passage_num_list:
                    cov_per_passage = sum([cov_per_pass_repl_dict[passage_num][repl_val]
                                           for repl_val in cov_per_pass_repl_dict[passage_num]])
                    list_of_passage_cov.append(cov_per_passage)
                    junc_per_pos_per_passage = sum([junc_per_pos_per_pass_repl_dict[passage_num][repl_val]
                                                    for repl_val in junc_per_pos_per_pass_repl_dict[passage_num]])
                    list_of_passage_junc_per_pos.append(
                        junc_per_pos_per_passage)

                list_of_passsage_sign = [
                    sign_per_pass_dict[passage_num] for passage_num in passage_num_list]
                list_of_passsage_sign_sum = [
                    sum_sign_per_pass_dict[passage_num] for passage_num in passage_num_list]
                list_of_passsage_num_points = [
                    num_points_per_pass_dict[passage_num] for passage_num in passage_num_list]
                out_file1.write(str(points_tuple[neighb.center_point_ind].x1) + '\t' +
                                str(points_tuple[neighb.center_point_ind].y1) + '\t' +
                                str(neighb.r[1]) + '\t' +
                                str(neighb.size) + '\t' +
                                str(neighb.pvalue) + '\t' +
                                '\t'.join(map(str, list_of_passsage_sign)) + '\n')
                out_file2.write(str(points_tuple[neighb.center_point_ind].x1) + '\t' +
                                str(points_tuple[neighb.center_point_ind].y1) + '\t' +
                                str(neighb.r[1]) + '\t' +
                                str(neighb.size) + '\t' +
                                str(neighb.pvalue) + '\t' +
                                '\t'.join(map(str, list_of_passsage_sign_sum)) + '\n')
                out_file3.write(str(points_tuple[neighb.center_point_ind].x1) + '\t' +
                                str(points_tuple[neighb.center_point_ind].y1) + '\t' +
                                str(neighb.r[1]) + '\t' +
                                str(neighb.size) + '\t' +
                                str(neighb.pvalue) + '\t' +
                                '\t'.join(map(str, list_of_passsage_num_points)) + '\n')
                # Cov and junct per pos
                out_file4.write(str(points_tuple[neighb.center_point_ind].x1) + '\t' +
                                str(points_tuple[neighb.center_point_ind].y1) + '\t' +
                                str(neighb.r[1]) + '\t' +
                                str(neighb.size) + '\t' +
                                str(neighb.pvalue) + '\t' +
                                '\t'.join(map(str, list_of_passage_cov + list_of_passage_junc_per_pos)) + '\n')

                # for replicate_num in replicate_num_list:
                #     out_file1.write(' ' + '\t' + \
                #             ' ' + '\t' + \
                #             ' '+ '\t' + \
                #             ' ' + '\t' + \
                #             ' ' + '\t' + \
                #             '\t'.join(map(str, [sign_per_pass_repl_dict[passage][replicate_num]
                #                                 for passage in sign_per_pass_repl_dict])) + '\n')
                # for replicate_num in replicate_num_list:
                #     out_file2.write(' ' + '\t' + \
                #             ' ' + '\t' + \
                #             ' '+ '\t' + \
                #             ' ' + '\t' + \
                #             ' ' + '\t' + \
                #             '\t'.join(map(str, [sum_sign_per_pass_repl_dict[passage][replicate_num]
                #                                 for passage in sum_sign_per_pass_repl_dict])) + '\n')
                # for replicate_num in replicate_num_list:
                #     out_file3.write(' ' + '\t' + \
                #             ' ' + '\t' + \
                #             ' '+ '\t' + \
                #             ' ' + '\t' + \
                #             ' ' + '\t' + \
                #             '\t'.join(map(str, [num_points_per_pass_repl_dict[passage][replicate_num]
                #                                 for passage in num_points_per_pass_repl_dict])) + '\n')
                # # Cov and junct per pos
                # for replicate_num in replicate_num_list:
                #     out_file4.write(' ' + '\t' + \
                #             ' ' + '\t' + \
                #             ' '+ '\t' + \
                #             ' ' + '\t' + \
                #             ' ' + '\t' + \
                #             '\t'.join(map(str, [cov_per_pass_repl_dict[passage][replicate_num]
                #                                 for passage in cov_per_pass_repl_dict])) + '\t' +\
                #             '\t'.join(map(str, [junc_per_pos_per_pass_repl_dict[passage][replicate_num]
                #                                 for passage in junc_per_pos_per_pass_repl_dict])) + '\n')
###
###

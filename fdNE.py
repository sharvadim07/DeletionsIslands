import logging
import argparse

parser = argparse.ArgumentParser(
    description='Fractal Dimension Neighbor Embedding (fdNE) algorithm in application to fitness of deletions in genome across groups (passages)')

# Paramaters for two dimension case with deletions
parser.add_argument('-d', '--del_tables', type=str,
                    help='Deletions table.', required=False)
parser.add_argument('-f', '--freq_filt', type=str,
                    help='Threshold for frequency filtering.', required=False)
parser.add_argument('-s', '--freq_filt_2', type=float,
                    help='Threshold 2 for frequency filtering.', required=False)
parser.add_argument('-l', '--del_length', type=int,
                    help='Deletions length threshold.', required=False)
parser.add_argument('-n', '--neighb_info', type=str,
                    help='File with neighborhood info.', required=False)
parser.add_argument('-c', '--only_coords', action='store_true',
                    help='Coords print flag.')
parser.add_argument('--do_two_step', action='store_true',
                    help='Do two step neighborhood identification?.')
parser.add_argument('--do_mod_step', action='store_true',
                    help='Do modified step neighborhood identification?.')

# Parameters for vectors of params from dynamics case
parser.add_argument('-p', '--par_tables', type=str,
                    help='Parameters tables.', required=False)
parser.add_argument('--phy_dyn', action='store_true',
                    help='Param table from phylo dynamics.')
parser.add_argument('--exp_dyn', action='store_true',
                    help='Param table from experimental dynamics.')
parser.add_argument('--single_cell', action='store_true',
                    help='Expression values table from single cell.')


# General params
parser.add_argument('-K', '--num_K_points', type=int,
                    help='Number of K random points.', required=False)  # 5
parser.add_argument('-A', '--A_value', type=str,
                    help='Value to rounding distances between all points and K points.', required=False)  # 10
# parser.add_argument('-Z', '--zone_size_percent', type=float,
#                     help='Zone size for min pvalue searching in percent from all space of points.', required=False)  # 100

def main():
    args = parser.parse_args()
    # 1 Type of input deletions table
    # Type of points is point on two dimensional plane
    if args.del_tables != None and args.freq_filt != None and args.del_length != None:
        import deleteions_table
        deleteions_table.deletions_table_processing(args, 'out_deletions_table/')
    # 2 Type of input params table
    # Type of points is vector of parameters on one MCMC iteration
    elif args.par_tables != None:
        import param_table
        if args.phy_dyn:
            param_table.param_table_processing(
                args, 'out_params_table/', 'phy_dyn')
        elif args.exp_dyn:
            param_table.param_table_processing(
                args, 'out_params_table/', 'exp_dyn')
        elif args.single_cell:
            param_table.param_table_processing(
                args, 'out_params_table/', 'single_cell')

if __name__ == "__main__":
    main()
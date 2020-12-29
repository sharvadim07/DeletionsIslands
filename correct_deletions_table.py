import argparse


parser = argparse.ArgumentParser(description='Correct deletions table program.')
parser.add_argument('-d','--del_tables', type=str,                   
                    help='Deletions table.', required=False)

args = parser.parse_args()

for file_name in args.del_tables.strip().split(','):
    import os
    with open(file_name, 'r') as in_file:
        with open('corr_' + os.path.basename(file_name), 'w') as out_file:
            for line in in_file:
                new_line_list = []
                line_list = line.strip().split('\t')
                for i, el in enumerate(line_list):
                    if i == 1:
                        el_list = el.strip().split('_')
                        new_line_list.append(el_list[0])
                        new_line_list.append(el_list[1])
                    else:
                        new_line_list.append(el)    
                out_file.write('\t'.join(new_line_list) + '\n')
                

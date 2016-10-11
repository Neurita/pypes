#!/usr/bin/env python
import os
import os.path as path
from multiprocessing import Pool

from docstamp.template import TextDocument
from docstamp.commands import call_command


template_file = 'rsn_multiple_regressions_template.m'

rsn_atlas_file = '/home/alexandre/data/std_brains/resting_state/allen/baseline/ALL_HC_unthresholded_tmaps_resampled.nii'
allen_rsn_idx = \
[['Basal_Ganglia_21', 21],
 ['Auditory_17',      17],
 ['Sensorimotor_7',    7],
 ['Sensorimotor_23',  23],
 ['Sensorimotor_24',  24],
 ['Sensorimotor_38',  38],
 ['Sensorimotor_56',  56],
 ['Sensorimotor_29',  29],
 ['Visual_46',        46],
 ['Visual_64',        64],
 ['Visual_67',        67],
 ['Visual_48',        48],
 ['Visual_39',        39],
 ['Visual_59',        59],
 ['Default_Mode_50',  50],
 ['Default_Mode_53',  53],
 ['Default_Mode_25',  25],
 ['Default_Mode_68',  68],
 ['Attentional_34',   34],
 ['Attentional_60',   60],
 ['Attentional_52',   52],
 ['Attentional_72',   72],
 ['Attentional_71',   71],
 ['Attentional_55',   55],
 ['Frontal_42',       42],
 ['Frontal_20',       20],
 ['Frontal_47',       47],
 ['Frontal_49',       49]]

output_dir     = '/home/alexandre/data/thomas/ica_out/8mm/fmri_no-grptemplate_noWMcor_30ICs/'
ica_param_file = '/home/alexandre/data/thomas/ica_out/8mm/fmri_no-grptemplate_noWMcor_30ICs/_ica_parameter_info.mat'
zscore_thr     = 2

n_cpus = 2


def run_icatb_multipleRegression(matlab_script_file):
    cur_path = path.realpath(os.curdir)

    os.chdir(output_dir)
    call_command('matlab', '-nodesktop -nojvm -nosplash'.split() +
                 ['-r', path.basename(matlab_script_file).split('.')[0]])
    os.chdir(cur_path)


if __name__ == '__main__':
    template_doc = TextDocument(template_file)

    template_args = {'output_dir': output_dir,
                     'ica_param_file': ica_param_file,
                     'rsn_atlas_file': rsn_atlas_file,
                     'zscore_threshold': zscore_thr,
                    }

    matlab_scripts = []

    # create the sets of arguments and render the matlab scripts from the template
    for name, idx in allen_rsn_idx:

        # add specific arguments to the list of arguments
        args = template_args.copy()
        args['output_basename']   = name
        args['template_ic_index'] = idx

        # the script path
        _ = template_doc.fill(doc_contents=args)
        matlab_script_file = path.join(output_dir,
                                       'icatb_multipleRegressions_to_{}.m'.format(name))

        # fill its arguments
        template_doc.render(matlab_script_file)

        # store its path to run it later
        matlab_scripts.append(matlab_script_file)

    # multiprocessing run
    pool = Pool(processes=n_cpus)

    result  = pool.map(run_icatb_multipleRegression, matlab_scripts)
    measures = [x for x in result if not x is None]
    #cleaned = np.asarray(cleaned)
    # not optimal but safe
    pool.close()
    pool.join()


    #
    # for name, idx in allen_rsn_idx:
    #     print('Processing {}.'.format(name))
    #     # if idx != 50:
    #     #     continue
    #     run_icatb_multipleRegression()
    #     contents['output_basename']   = name
    #     contents['template_ic_index'] = idx
    #
    #     _ = template_doc.fill(doc_contents=contents)
    #     matlab_script_file = path.join(output_dir,
    #                                    'icatb_multipleRegressions_to_{}.m'.format(name))
    #     template_doc.render(matlab_script_file)
    #
    #     os.chdir(output_dir)
    #     call_command('matlab', '-nodesktop -nojvm -nosplash'.split() +
    #                  ['-r', path.basename(matlab_script_file).split('.')[0]])
    #     os.chdir(cur_path)

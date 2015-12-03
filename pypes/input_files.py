"""
Workflows to grab input file structures.
"""
import os
import os.path as op

import nipype.pipeline.engine as pe
from   nipype.interfaces.utility import IdentityInterface
from   nipype.interfaces.io import SelectFiles

from ._utils import _check_list


def subject_session_input(base_dir, session_names, file_names, subject_ids=None,
                          wf_name=None):
    """ A workflow of IdentityInterface->SelectFiles for the case where
    you have a {subject_id}/{session_id}/{image_file} dataset structure.

    Parameters
    ----------
    base_dir: str
        Path to the working directory of the workflow

    session_names: list of str
        Example: ['session_0']
        If None or '', will remove the 'session_id' option from the information source.

    file_names: Dict[str -> str]
        A dictionary that relates the `select` node keynames and the
        file name.
        Example: {'anat': 'anat_hc.nii.gz',       'pet': 'pet_fdg.nii.gz'},
        Example: {'anat': 'anat_1/mprage.nii.gz', 'rest': 'rest_1/rest.nii.gz'},

    subject_ids: list of str
        Use this if you want to limit the analysis to certain subject IDs.
        If `None` will pick the folders from os.listdir(data_dir).

    wf_name: str
        Name of the workflow
        Default: "subject_session_files"

    Nipype Outputs
    --------------
    select.{file_name_without_extension}: path to existing file
        Will give you the path of the {file_name}s.

    infosrc.subject_id: str
        Will give you the 'subject_id's.

    infosrc.session_id: str
        Will give you the 'session_id's.


    Returns
    -------
    wf: nipype Workflow
    """
    # check default values
    if wf_name is None:
        wf_name = "subject_session_files"

    # Check the subject ids
    subj_ids = _check_list(subject_ids)
    if subj_ids is None:
        subj_ids = [op.basename(p) for p in os.listdir(base_dir)]

    # the fields and its values for the fields_iterables
    fields  = [('subject_id', subj_ids)]
    dirpath = '{{subject_id}}'
    if session_names is not None:
        fields.append(('session_id', session_names))
        dirpath = op.join(dirpath, '{{session_id}}')

    files = {'{}'.format(f): op.join(dirpath, '{}').format(file_names[f]) for f in file_names}

    return input_file_wf(work_dir=base_dir,
                         data_dir=base_dir,
                         field_iterables=fields,
                         file_templates=files,
                         wf_name=wf_name)


def input_file_wf(work_dir, data_dir, field_iterables, file_templates, wf_name="input_files"):
    """ A workflow of IdentityInterface->SelectFiles for the case where
    you have a {subject_id}/{session_id}/{image_file} dataset structure.

    Parameters
    ----------
    work_dir: str
        Path to the working directory of the workflow

    data_dir: str

    field_iterables: List of 2-tuples (str, iterable)
        Example: [('session_id', session_names),
                  ('subject_id', subject_ids),]
        This will be input to an IdentityInterface

    file_templates: Dict[str -> str]
        Example: {'anat': '{subject_id}/{session_id}/anat_hc.nii.gz',
                  'pet': '{subject_id}/{session_id}/pet_fdg.nii.gz',
                 }

    wf_name: str
        Name of the workflow

    Nipype Outputs
    --------------
    select.{template_key}: path to existing file
        Will give you the path of the {file_name}s taken from `file_templates`.

    infosrc.{field_name}: str
        Will give you the value of the field from `field_iterables`.

    Returns
    -------
    wf: nipype Workflow
    """
    # Input workflow
    wf = pe.Workflow(name=wf_name, base_dir=work_dir)

    # Infosource - a function free node to iterate over the list of subject names
    field_names = [field[0] for field in field_iterables]

    infosource = pe.Node(IdentityInterface(fields=field_names), name="infosrc")
    infosource.iterables = field_iterables

    # SelectFiles
    select = pe.Node(SelectFiles(file_templates, base_directory=data_dir), name="select")

    # Connect, e.g., 'infosrc.subject_id' to 'select.subject_id'
    wf.connect([(infosource, select, [(field, field) for field in field_names])])

    return wf


def get_input_file_name(input_node, fname_key):
    """ Return the name of the file given by the node key `fname_key` in `input_node`.

    Parameters
    ----------
    input_node: nipype Node
        a node with a file input interface (SelectFiles, for now).

    fname_key: str
        The key that is used to access the file path using `input_node`.
        Example: 'anat'

    Returns
    -------
    filename: str
        The base input file path from the input node.
    """
    if isinstance(input_node.interface, SelectFiles):
        try:
            fname = input_node.interface._templates[fname_key]
        except:
            raise AttributeError("Could not find a SelectFiles node called 'select' in main workflow.")
        else:
            return fname
    else:
        raise NotImplementedError('`get_input_file_name` has not been implemented for nodes'
                                  ' of type {}.'.format(type(input_node.interface)))

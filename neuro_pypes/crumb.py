# -*- coding: utf-8 -*-
"""
A nipype.SelectFiles node based on hansel.Crumb
"""
import os.path as op
from warnings import warn

from hansel import Crumb

from nipype.interfaces.base import (traits,
                                    DynamicTraitedSpec,
                                    Undefined, BaseInterfaceInputSpec)
from nipype.interfaces.io import IOBase, add_traits
from nipype.utils.filemanip import list_to_filename
from nipype.utils.misc import human_order_sorted

from ._utils import get_values_map_keys


class DataCrumbInputSpec(DynamicTraitedSpec, BaseInterfaceInputSpec):
    sort_filelist = traits.Bool(True, usedefault=True,
                                desc='Sort the filelist that matches the template. Crumb always sort its outputs.')
    raise_on_empty = traits.Bool(True, usedefault=True,
                                 desc="Raise an exception if a template pattern matches no files.")
    force_lists = traits.Either(traits.Bool(), traits.List(traits.Str()),
                                default=False, usedefault=True,
                                desc=("Whether to return outputs as a list even when only one file "
                                      "matches the template. Either a boolean that applies to all "
                                      "output fields or a list of output field names to coerce to "
                                      " a list"))


class DataCrumb(IOBase):
    """ Use Crumb from hansel to select input files."""
    input_spec  = DataCrumbInputSpec
    output_spec = DynamicTraitedSpec
    _always_run = True

    def __init__(self, crumb, templates, **kwargs):
        """Create an instance with specific input fields.
        Parameters
        ----------
        crumb: hansel.Crumb
            If you are using a relative crumb path use a first argument as
            base directory. This argument will be exposed as an input.
            Example:
            {base_dir}/data/raw/{subj_id}...

        templates : dict[str] -> list of 2-tuples
            Mapping from string keys to list of crumb arguments in crumb_path
            that must be replaced to complete the file crumb path.

            The keys become output fields on the interface.

            At runtime, the values of the interface inputs will be
            plugged into these templates, and the resulting strings will be
            used to select files.
        """
        super(DataCrumb, self).__init__(**kwargs)

        # Infer the infields and outfields from the template
        if not Crumb.is_valid(crumb):
            raise ValueError('Crumb {} is not valid.'.format(crumb))

        self._crumb = crumb

        files_args = get_values_map_keys(templates)
        undef_args = [name for name in list(crumb.all_args()) if name not in files_args]
        self._infields  = undef_args

        self._outfields = []
        self._templates = []
        if templates:
            self._outfields = list(templates)
            self._templates = templates

        # Add the dynamic input fields
        undefined_traits = {}
        for field in self._infields:
            self.inputs.add_trait(field, traits.Any)
            undefined_traits[field] = Undefined
        self.inputs.trait_set(trait_change_notify=False, **undefined_traits)

    def _add_output_traits(self, base):
        """Add the dynamic output fields"""
        out_fields = list(self._templates.keys()) + list(self._infields)
        return add_traits(base, out_fields)

    def _list_outputs(self):
        """Find the files and expose them as interface outputs."""
        outputs = {}
        info = dict([(k, v) for k, v in list(self.inputs.__dict__.items())
                     if k in self._infields])

        # check if the crumb is not absolute or if in info we have the parameter for the base directory
        if not self._crumb.isabs():
            first_arg_name, _ = self._crumb._first_open_arg()
            if first_arg_name not in info:
                raise KeyError('Crumb path is not absolute and could not find input for {}.'.format(first_arg_name))
            elif not op.isabs(info[first_arg_name]):
                raise IOError('Expected an absolute path for {} argument in {} but got {}.'.format(first_arg_name,
                                                                                                   self._crumb,
                                                                                                   info[first_arg_name],
                                                                                                   ))
        force_lists = self.inputs.force_lists
        if isinstance(force_lists, bool):
            force_lists = self._outfields if force_lists else []
        bad_fields = set(force_lists) - set(self._outfields)
        if bad_fields:
            bad_fields = ", ".join(list(bad_fields))
            plural = "s" if len(bad_fields) > 1 else ""
            verb = "were" if len(bad_fields) > 1 else "was"
            msg = ("The field%s '%s' %s set in 'force_lists' and not in "
                   "'templates'.") % (plural, bad_fields, verb)
            raise ValueError(msg)

        # loop over the crumb arguments to fill self_crumb
        crumb_info = {k: v for k, v in info.items() if k in self._crumb.open_args()}
        ocrumb = self._crumb.replace(**crumb_info)

        # check again if crumb path is absolute
        if not ocrumb.isabs():
            raise ValueError('Expected a Crumb with an absolute path, got {}.'.format(ocrumb))

        if not ocrumb.exists():
            raise IOError('Expected an existing Crumb path, got {}.'.format(ocrumb))

        # loop over all the ouput items and fill them with the info in templates
        for field, template in self._templates.items():

            # Fill in the template and glob for files
            focrumb  = ocrumb.replace(**dict(template))

            if list(focrumb.open_args()):
                raise ValueError('Expected a full specification of the Crumb path by now, got {}.'.format(focrumb))

            filelist = [cr.path for cr in focrumb.unfold()]
            # Handle the case where nothing matched
            if not filelist:
                msg = "No files were found unfolding %s crumb path: %s" % (
                    field, focrumb)
                if self.inputs.raise_on_empty:
                    raise IOError(msg)
                else:
                    warn(msg)

            # Possibly sort the list
            if self.inputs.sort_filelist:
                filelist = human_order_sorted(filelist)

            # Handle whether this must be a list or not
            if field not in force_lists:
                filelist = list_to_filename(filelist)

            outputs[field] = filelist

            # add the crumb argument values for output
            for arg_name in focrumb.all_args():
                outputs[arg_name] = focrumb[arg_name][0]

        return outputs

# -*- coding: utf-8 -*-
"""
Utilities for the CLI functions.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os.path as path
import re
from typing import Iterable, Tuple

import click

import hansel

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
UNKNOWN_OPTIONS = dict(allow_extra_args=True,
                       ignore_unknown_options=True)

ExistingDirPath = click.Path(exists=True, file_okay=False, resolve_path=True)
ExistingFilePath = click.Path(exists=True, dir_okay=False, resolve_path=True)
UnexistingFilePath = click.Path(dir_okay=False, resolve_path=True)


def check_not_none(ctx, param, value):
    if value is None:
        raise click.BadParameter('got {}.'.format(value))
    return value


class RegularExpression(click.ParamType):
    name = 'regex'

    def convert(self, value, param, ctx):
        try:
            rex = re.compile(value, re.IGNORECASE)
        except ValueError:
            self.fail('%s is not a valid regular expression.' % value, param, ctx)
        else:
            return rex


class CrumbPath(click.ParamType):
    name = 'crumb'

    def convert(self, value, param, ctx):
        try:
            cr = hansel.Crumb(path.expanduser(value), ignore_list=['.*'])
        except ValueError:
            self.fail('%s is not a valid crumb path.' % value, param, ctx)
        else:
            return cr


def echo_list(alist):
    for i in alist:
        click.echo(i)


def _print_values_map_as_csv(list_of_lists):
    for values in list_of_lists:
        click.echo(','.join([item[1] for item in values]))


def _get_file_pairs(background: hansel.Crumb, foreground: hansel.Crumb) -> Iterable[Tuple[str, [str, None]]]:
    if not foreground:
        return ((bg, None) for bg in background.ls(make_crumbs=False))

    if not background.has_crumbs() and not foreground.has_crumbs():
        return [(background.ls(make_crumbs=False)[0], foreground.ls(make_crumbs=False)[0])]

    if foreground.has_crumbs() and background.has_crumbs():
        background_args = tuple(background.open_args())
        foreground_args = tuple(foreground.open_args())
        if background_args != foreground_args:
            raise click.BadParameter('Make sure that background and foreground have the same crumb arguments. '
                                     'Got {}, and {}.'.format(background_args, foreground_args))

        try:
            intersect_values = hansel.intersection(background, foreground)
        except KeyError as ke:
            raise click.BadParameter('Make sure that background and foreground have the same crumb arguments. '
                                     'Got {}, and {}.'.format(background_args, foreground_args)) from ke
        else:
            return (
                (background.build_paths(crumb_args, make_crumbs=False)[0],
                 foreground.build_paths(crumb_args, make_crumbs=False)[0]) for crumb_args in intersect_values)

    bg_files = []
    fg_files = []
    if background.has_crumbs():
        bg_files = list(background.ls(make_crumbs=False))
        fg_files = [foreground.ls(make_crumbs=False)[0]] * len(bg_files)

    if not foreground.has_crumbs():
        fg_files = list(foreground.ls(make_crumbs=False))
        bg_files = [background.ls(make_crumbs=False)[0]] * len(fg_files)

    return zip(bg_files, fg_files)

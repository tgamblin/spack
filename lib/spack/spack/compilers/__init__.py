##############################################################################
# Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# This file is part of Spack.
# Created by Todd Gamblin, tgamblin@llnl.gov, All rights reserved.
# LLNL-CODE-647188
#
# For details, see https://github.com/llnl/spack
# Please also see the LICENSE file for our notice and the LGPL.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License (as
# published by the Free Software Foundation) version 2.1, February 1999.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
# conditions of the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##############################################################################
"""This module contains functions related to finding compilers on the
system and configuring Spack to use multiple compilers.
"""
import imp
import platform

from llnl.util.lang import list_modules
from llnl.util.filesystem import join_path

import spack
import spack.error
import spack.spec
import spack.config
import spack.architecture

from spack.util.naming import mod_to_class

_imported_compilers_module = 'spack.compilers'
_path_instance_vars = ['cc', 'cxx', 'f77', 'fc']
_other_instance_vars = ['modules', 'operating_system']
_cache_config_file = []

# TODO: customize order in config file
if platform.system() == 'Darwin':
    _default_order = ['clang', 'gcc', 'intel']
else:
    _default_order = ['gcc', 'intel', 'pgi', 'clang', 'xlc', 'nag']


def _auto_compiler_spec(function):
    def converter(cspec_like, *args, **kwargs):
        if not isinstance(cspec_like, spack.spec.CompilerSpec):
            cspec_like = spack.spec.CompilerSpec(cspec_like)
        return function(cspec_like, *args, **kwargs)
    return converter


def _to_dict(compiler):
    """Return a dict version of compiler suitable to insert in YAML."""
    d = {}
    d['spec'] = str(compiler.spec)
    d['paths'] = dict((attr, getattr(compiler, attr, None))
                      for attr in _path_instance_vars)
    d['flags'] = dict((fname, fvals) for fname, fvals in compiler.flags)
    d['operating_system'] = str(compiler.operating_system)
    d['modules'] = compiler.modules if compiler.modules else []

    if compiler.alias:
        d['alias'] = compiler.alias

    return {'compiler': d}


def get_compiler_config(scope=None, init_config=True):
    """Return the compiler configuration for the specified architecture.
    """
    def init_compiler_config():
        """Compiler search used when Spack has no compilers."""
        compilers = find_compilers()
        compilers_dict = []
        for compiler in compilers:
            compilers_dict.append(_to_dict(compiler))
        spack.config.update_config('compilers', compilers_dict, scope=scope)

    config = spack.config.get_config('compilers', scope=scope)
    # Update the configuration if there are currently no compilers
    # configured.  Avoid updating automatically if there ARE site
    # compilers configured but no user ones.
    if not config and init_config:
        if scope is None:
            # We know no compilers were configured in any scope.
            init_compiler_config()
            config = spack.config.get_config('compilers', scope=scope)
        elif scope == 'user':
            # Check the site config and update the user config if
            # nothing is configured at the site level.
            site_config = spack.config.get_config('compilers', scope='site')
            if not site_config:
                init_compiler_config()
                config = spack.config.get_config('compilers', scope=scope)
        return config
    elif config:
        return config
    else:
        return []  # Return empty list which we will later append to.


def add_compilers_to_config(compilers, scope=None, init_config=True):
    """Add compilers to the config for the specified architecture.

    Arguments:
      - compilers: a list of Compiler objects.
      - scope:     configuration scope to modify.
    """
    compiler_config = get_compiler_config(scope, init_config)
    for compiler in compilers:
        compiler_config.append(_to_dict(compiler))
    global _cache_config_file
    _cache_config_file = compiler_config
    spack.config.update_config('compilers', compiler_config, scope)


@_auto_compiler_spec
def remove_compiler_from_config(compiler_spec, scope=None):
    """Remove compilers from the config, by spec.

    Arguments:
      - compiler_specs: a list of CompilerSpec objects.
      - scope:          configuration scope to modify.
    """
    # Need a better way for this
    global _cache_config_file

    compiler_config = get_compiler_config(scope)
    config_length = len(compiler_config)

    filtered_compiler_config = [
        comp for comp in compiler_config
        if spack.spec.CompilerSpec(comp['compiler']['spec']) != compiler_spec]

    # Update the cache for changes
    _cache_config_file = filtered_compiler_config
    if len(filtered_compiler_config) == config_length:  # No items removed
        CompilerSpecInsufficientlySpecificError(compiler_spec)
    spack.config.update_config('compilers', filtered_compiler_config, scope)


def all_compilers_config(scope=None, init_config=True):
    """Return a set of specs for all the compiler versions currently
       available to build with.  These are instances of CompilerSpec.
    """
    # Get compilers for this architecture.
    # Create a cache of the config file so we don't load all the time.
    global _cache_config_file
    if not _cache_config_file:
        _cache_config_file = get_compiler_config(scope, init_config)
        return _cache_config_file
    else:
        return _cache_config_file


def all_compilers(scope=None, init_config=True):
    """Return compiler specs from the merged config."""
    cspecs = []
    for s in all_compilers_config(scope, init_config):
        cspec = spack.spec.CompilerSpec(s['compiler']['spec'])

        # ensure versions read from config are exact
        cspec.versions = cspec.versions.exactly
        cspecs.append(cspec)
    return cspecs


def default_compiler():
    versions = []
    for name in _default_order:
        versions = find(name)
        if versions:
            break
    else:
        raise NoCompilersError()

    return sorted(versions)[-1]


def find_compilers(*paths):
    """Return a list of compilers found in the suppied paths.
       This invokes the find_compilers() method for each operating
       system associated with the host platform, and appends
       the compilers detected to a list.
    """
    # Find compilers for each operating system class
    oss = all_os_classes()
    compiler_lists = []
    for o in oss:
        compiler_lists.extend(o.find_compilers(*paths))
    return compiler_lists


def supported_compilers():
    """Return a set of names of compilers supported by Spack.

       See available_compilers() to get a list of all the available
       versions of supported compilers.
    """
    return sorted(name for name in list_modules(spack.compilers_path))


@_auto_compiler_spec
def supported(compiler_spec):
    """Test if a particular compiler is supported."""
    return compiler_spec.name in supported_compilers()


@_auto_compiler_spec
def find(compiler_spec, scope=None):
    """Return specs of available compilers that match the supplied
       compiler spec.  Return an empty list if nothing found."""
    return [c for c in all_compilers(scope) if c.satisfies(compiler_spec)]


@_auto_compiler_spec
def compilers_for_spec(compiler_spec, scope=None, **kwargs):
    """This gets all compilers that satisfy the supplied CompilerSpec.
       Returns an empty list if none are found.
    """
    platform = kwargs.get('platform', None)
    config = all_compilers_config(scope)

    def get_compilers(cspec):
        compilers = []

        for items in config:
            if items['compiler']['spec'] != str(cspec):
                continue
            items = items['compiler']

            if not ('paths' in items and
                    all(n in items['paths'] for n in _path_instance_vars)):
                raise InvalidCompilerConfigurationError(cspec)

            cls  = class_for_compiler_name(cspec.name)

            compiler_paths = []
            for c in _path_instance_vars:
                compiler_path = items['paths'][c]
                if compiler_path != 'None':
                    compiler_paths.append(compiler_path)
                else:
                    compiler_paths.append(None)

            mods = items.get('modules')
            if mods == 'None':
                mods = []

            os = None
            if 'operating_system' in items:
                os = spack.architecture._operating_system_from_dict(
                    items['operating_system'], platform)

            alias = items.get('alias', None)

            compiler_flags = items.get('flags', {})

            compilers.append(
                cls(cspec, os, compiler_paths, mods, alias, **compiler_flags))

        return compilers

    matches = set(find(compiler_spec, scope))
    compilers = []
    for cspec in matches:
        compilers.extend(get_compilers(cspec))
    return compilers


@_auto_compiler_spec
def compiler_for_spec(compiler_spec, arch):
    """Get the compiler that satisfies compiler_spec.  compiler_spec must
       be concrete."""
    operating_system = arch.platform_os
    assert(compiler_spec.concrete)

    compilers = [
        c for c in compilers_for_spec(compiler_spec, platform=arch.platform)
        if c.operating_system == operating_system]
    if len(compilers) < 1:
        raise NoCompilerForSpecError(compiler_spec, operating_system)
    if len(compilers) > 1:
        raise CompilerSpecInsufficientlySpecificError(compiler_spec)
    return compilers[0]


def class_for_compiler_name(compiler_name):
    """Given a compiler module name, get the corresponding Compiler class."""
    assert(supported(compiler_name))

    file_path = join_path(spack.compilers_path, compiler_name + ".py")
    compiler_mod = imp.load_source(_imported_compilers_module, file_path)
    cls = getattr(compiler_mod, mod_to_class(compiler_name))

    # make a note of the name in the module so we can get to it easily.
    cls.name = compiler_name

    return cls


def all_os_classes():
    """
    Return the list of classes for all operating systems available on
    this platform
    """
    classes = []

    platform = spack.architecture.platform()
    for os_class in platform.operating_sys.values():
        classes.append(os_class)

    return classes


def all_compiler_types():
    return [class_for_compiler_name(c) for c in supported_compilers()]


class InvalidCompilerConfigurationError(spack.error.SpackError):

    def __init__(self, compiler_spec):
        super(InvalidCompilerConfigurationError, self).__init__(
            "Invalid configuration for [compiler \"%s\"]: " % compiler_spec,
            "Compiler configuration must contain entries for all compilers: %s"
            % _path_instance_vars)


class NoCompilersError(spack.error.SpackError):
    def __init__(self):
        super(NoCompilersError, self).__init__(
            "Spack could not find any compilers!")


class NoCompilerForSpecError(spack.error.SpackError):
    def __init__(self, compiler_spec, target):
        super(NoCompilerForSpecError, self).__init__(
            "No compilers for operating system %s satisfy spec %s"
            % (target, compiler_spec))


class CompilerSpecInsufficientlySpecificError(spack.error.SpackError):
    def __init__(self, compiler_spec):
        super(CompilerSpecInsufficientlySpecificError, self).__init__(
            "Multiple compilers satisfy spec %s" % compiler_spec)

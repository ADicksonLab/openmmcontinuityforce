from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
contforceplugin_header_dir = '@CONTFORCEPLUGIN_HEADER_DIR@'
contforceplugin_library_dir = '@CONTFORCEPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = ['-std=c++11']
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_contforceplugin',
                      sources=['ContForcePluginWrapper.cpp'],
                      libraries=['OpenMM', 'ContForcePlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), contforceplugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), contforceplugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='contforceplugin',
      version='1.0',
      py_modules=['contforceplugin'],
      ext_modules=[extension],
     )

from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
seekr2plugin_header_dir = '@SEEKR2PLUGIN_HEADER_DIR@'
seekr2plugin_library_dir = '@SEEKR2PLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = ['-std=c++11']
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_seekr2plugin',
                      sources=['Seekr2PluginWrapper.cpp'],
                      libraries=['OpenMM', 'Seekr2Plugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), seekr2plugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), seekr2plugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='seekr2plugin',
      version='0.1',
      py_modules=['seekr2plugin'],
      ext_modules=[extension],
     )

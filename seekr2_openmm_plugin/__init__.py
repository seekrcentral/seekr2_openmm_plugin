"""
seekr2_openmm_plugin
An MMVT OpenMM plugin for SEEKR2.
"""

# Add imports here
from .mmvtplugin import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

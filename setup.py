"""
MicArTools
A toolset for navigating and analyzing microarray datasets
alias: micartools

Copyright (C) 2019  Jordan A. Berg
jordan <dot> berg <at> biochem <dot> utah <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from setuptools import setup
import re

with open('micartools/__init__.py', 'r') as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        fd.read(), re.MULTILINE).group(1)

setup(
    name = 'MicArTools',
    version = version,
    description = 'A toolset for navigating and analyzing microarray datasets',
    long_description = open('README.md').read(),
    author = 'Jordan Berg',
    author_email = 'jordan.berg@biochem.utah.edu',
    url = 'https://github.com/j-berg/micartools',
    packages = ['micartools'],
    package_dir = {'micartools': 'micartools'},
    license = 'GPL-3.0',
    zip_safe = False,
    install_requires=[
          'GEOparse',
          'pandas',
          'scipy',
          'seaborns',
          'plotly',
          'matplotlib',
          'sklearn',
          'mpl_toolkits',
          
      ],

    classifiers=[
        'Development Status :: Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
)

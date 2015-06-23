#!/usr/bin/env python

from distutils.core import setup

setup(name='OptiType',
      version='1.0',
      description='Precision HLA typing from next-generation sequencing data',
      author='Andras Szolek, Benjamin Schubert, Christopher Mohr',
      url='https://github.com/FRED-2/OptiType',
      license='LICENSE',
      # provided packages
      packages=['optitype'],
      # data files required by optitype
      package_data={'optitype': ['config.ini', 'data/*']},
      # provided binary packages
      scripts=['bin/optitype'],
      # packages that we depend on
      install_requires=[
        'Coopr>=3.5,<4.0',
        'argparse==1.2.1',
        'biopython==1.65',
        'matplotlib==1.4.3',
        'wsgiref==0.1.2',
        'pandas==0.16.2',
        'tables==3.2.0',
      ]
      )

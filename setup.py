#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name='meneval',
      url='',
      license='',
      description='Validation of Meneco reactions for gapfilling',
      long_description='',
      author='Pauline Hamon-Giraud',
      author_email='pauline.hamon-giraud@irisa.fr',
      python_requires='>=3.6',
      packages=['meneval'],
      install_requires=['meneco', 'biopython', 'padmet'],
      entry_points={
          'console_scripts': [
              'meneval = meneval.__main__:main'
          ]
      },
      )

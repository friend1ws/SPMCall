#!/usr/bin/env python

from distutils.core import setup

setup(name='mutrans',
      version='0.1.0',
      description='Python tools for extracting splicing-causing mutations',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/Genomon-Project/mutrans.git',
      package_dir = {'': 'lib'},
      packages=['mutrans'],
      scripts=['mutrans'],
      license='GPL-3'
     )


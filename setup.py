from setuptools import setup

__author__ = 'Stephen G. Gaffney'

setup(name='cnv_pipeline',
      version='0.1',
      description='Run ADTEx and saasCNV.',
      url='http://github.com/sggaffney',
      author='Stephen G. Gaffney',
      author_email='stephen.gaffney@yale.edu',
      license='GPLv3',
      packages=['cnv_pipeline'],
      package_data={'cnv_pipeline': ['config.ini', '*.R']},
      install_requires=[
          'pandas', 'feather-format==0.2.0', 'configparser'
      ],
      entry_points={'console_scripts': ['run_cnv = cnv_pipeline.pipeline:main']},
      zip_safe=False,
      )

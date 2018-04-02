from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='gispy',
      version='0.0.4',
      description='GIS functions with GDAL/OGR',
      long_description=readme(),
      classifiers=[
            'Development Status :: 3 - Alpha',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Natural Language :: English',
            'Topic :: Scientific/Engineering :: GIS',
      ],
      keywords='GIS, raster, vector',
      url='https://github.com/konradhafen/gispy',
      author='Konrad Hafen',
      author_email='khafen74@gmail.com',
      license='GPLv3',
      packages=['gispy'],
      install_requires=[
            'gdal',
            'numpy',
            'scipy',
      ],
      include_package_data=True,
      zip_safe=False)
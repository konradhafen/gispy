from setuptools import setup

setup(name='gispy',
      version='0.0.2',
      description='GIS functions with GDAL/OGR',
      url="https://github.com/konradhafen/gispy",
      author='Konrad Hafen',
      author_email='khafen74@gmail.com',
      license='GPLv3',
      packages=['gispy'],
      install_requires=[
            'gdal',
            'numpy',
            'scipy',
            'rasterstats',
      ],
      zip_safe=False)
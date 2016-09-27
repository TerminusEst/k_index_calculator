from setuptools import setup

def readme():
	with open("README.rst") as f:
		return f.read()

setup(name='k_index_calculator',
      version='0.2.3',
      description='Python module which calculates the K-Index of a geomagnetic time series.',
      long_description='Calculates the K-Index of geomagnetic time series using the FMI method',
      keywords='geomagnetic k-index space weather',
      url='https://github.com/TerminusEst/k_index_calculator',
      author='Sean Blake',
      author_email='blakese@tcd.ie',
      license='MIT',
      packages=['k_index_calculator'],
      install_requires=[
            'datetime', 'numpy', 'scipy', 'matplotlib'
            ],
      include_package_data=True,
      zip_safe=False)


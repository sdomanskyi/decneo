"""
setup.py for decneo package
"""
from setuptools import setup, find_packages
from codecs import open
from os import path
here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description=f.read()

setup(
    name='decneo',
    packages=find_packages(),
    version='1.0.6',
    description='Comberons from single cell transcriptomics in endothelial cells',
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    author='S. Domanskyi, A. Hakansson, M. Meng, J. S. Graff Zivin, C. Piermarocchi, G. Paternostro, N. Ferrara',
    author_email='s.domanskyi@gmail.com',
    license='MIT License',
    url='https://github.com/sdomanskyi/decneo',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: Unix',
        'Topic :: Education',
        'Topic :: Utilities',
        ],
    python_requires='>=3',
    install_requires=[
        'numpy>=1.19.1',
        'pandas>=1.0.1',
        'patsy>=0.5.1',
        'xlrd>=1.2.0',
        'openpyxl>=3.0.3',
        'tables>=3.6.1',
        'scipy>=1.4.1',
        'matplotlib>=3.1.3',
        'scikit-learn>=0.22.1',
        'networkx>=2.4',
        'adjustText>=0.7.3'],
    zip_safe=False
)

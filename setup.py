import pathlib
from setuptools import setup

VERSION = '1.3.0'
DESCRIPTION = 'An alignment-independent tool for the study of pathogen sequence diversity at any given rank of taxonomy lineage'

HERE = pathlib.Path(__file__).parent
README = (HERE/ "README.md").read_text()

setup(
    name='uniqmin',
    version=VERSION,
    description=DESCRIPTION,
    long_description=README,
    long_description_content_type="text/markdown",
    author="Chong LC",
    author_email="<lichuinchong@gmail.com>",
    keywords=("virus, pathogen, sequence diversity, alignment-independent, minimal set, proteome"),
    project_urls = 
      {"Github": "https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin"},
    packages=["uniqmin"],
    install_requires=['biopython==1.79', 'pandas>=1.3.5', 'pyahocorasick==1.4.1'],
    entry_points={
        "console_scripts": [
            "uniqmin=uniqmin.__main__:main",
        ]
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: Unix",
        "License :: OSI Approved :: MIT License", 
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics"        
    ]
)


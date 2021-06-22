#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import find_packages
from setuptools import setup

short_descr = "Segmentation and Tracking algorithm from Contact-dependent cell communications drive morphological invariance during ascidian embryogenesis."
readme = open('README.md').read()

# find packages
pkgs = find_packages('src')

setup_kwds = dict(
    name='astec',
    version="0.0.1",
    description=short_descr,
    long_description=readme,
    author="Gregoire Malandain",
    author_email="gregoire.malandain@inria.fr",
    url='https://github.com/astec-segmentation/astec',
    license='GPL',
    zip_safe=False,

    packages=pkgs,
    package_dir={'': 'src'},
    python_requires='>=3.7',
    setup_requires=[],
    install_requires=[],
    tests_require=[],
    entry_points={
        'console_scripts': [
            'astec_fusion=astec.bin.astec_fuse:main',
            'astec_mars=astec.bin.astec_mars:main',
            'astec_manual_correction=astec.bin.astec_manualcorrection:main',
            'astec_astec=astec.bin.astec_astec:main',
            'astec_postcorrection=astec.bin.astec_postcorrection:main',
            'astec_intraregistration=astec.bin.astec_intraregistration:main',
            'astec_embryoproperties=astec.bin.astec_embryoproperties:main',
            'astec_naming=astec.bin.astec_naming:main',
        ],
    },
    keywords='',

    test_suite='pytest',
    )

setup(**setup_kwds)

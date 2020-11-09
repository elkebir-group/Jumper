#! /usr/bin/env python
# -*- coding: utf-8 -*-
import setuptools

with open('README.md') as f:
    long_description = f.read()

setuptools.setup(
    name='jumper',
    packages=["jumper"],
    description="Discontinuous transcript assembly for coronaviruses",
    long_description=long_description,
    long_description_content_type='text/markdown',
    version='0.1.0',
    url='http://github.com/elkebir-group/jumper',
    author='Palash Sashittal and Chuanyi Zhang',
    author_email='sashitt2@illinois.edu',
    python_requires='>=3.6',
    # entry_points={'console_scripts': [
    #     'jumper = jumper.jumper:main_cli',
    #     'jumper_simulate = jumper.jumper_simulate:main_cli',
    # ]},
    scripts=[
        'scripts/jumper',
        'scripts/jumper_simulate',
    ],
    install_requires=[
        "pysam",
        "pandas",
        "gurobi",
    ],
)


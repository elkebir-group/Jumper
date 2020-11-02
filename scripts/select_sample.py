#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Thu June 18 2020
# @author: Chuanyi Zhang

"""
Select samples with depth profile close to canonical junction spanning
"""

import numpy as np
import pandas as pd
import math as m
import argparse
import sys
from os import listdir
from os.path import isfile, join


def max_pool(arr, size):
    steps = m.ceil(len(arr) / size)
    return np.array([np.max(arr[i*size:(i+1)*size]) for i in range(steps)])


def get_rss(arr):
    arr = arr / arr.max()
    # fit an exponential in pooled array
    pooled = max_pool(arr, 1000)
    x = np.array(range(15, len(pooled)))
    y = pooled[15:]
    # add 0.1 to prevent log(0)
    coeff = np.polyfit(x, np.log10(y+0.1), 1, w=np.sqrt(y+0.1))
    # calculate residual sum of square in original array
    coeff[0] /= 1000
    f = np.poly1d(coeff)
    return np.square(arr[15000:] - np.power(10, f(range(15000, len(arr))))).sum(), coeff


def by_coverage(arr):
    return np.mean(arr) >= 50 and np.mean(arr >= 20) >= 0.8


def main(args):
    depth = dict()
    path = args.folder
    for file in listdir(path):
        if file.endswith("depth"):
            sample = file.split('/')[-1].split('.')[0]
            depth[sample] = pd.read_csv(join(path, file),
                                        delimiter='\t',
                                        names=['chrom', 'pos', 'depth'])["depth"]

    candidate_samples = list(depth.keys())

    if args.cov:
        candidate_samples = [sample for sample in candidate_samples if by_coverage(depth[sample])]

    if args.fit_step:
        rss = dict()
        poly = dict()
        for sample in candidate_samples:
            dep = depth[sample]
            if len(dep) == 0:
                continue
            try:
                err, coeff = get_rss(dep)
                rss[sample] = err
                poly[sample] = coeff
            except ValueError:
                print(sample)
                continue

        selected_samples = dict()
        for sample in rss.keys():
            if rss[sample] <= 3500 and 0.000015 <= poly[sample][0]:
                selected_samples[sample] = (rss[sample], poly[sample])
        candidate_samples = list(selected_samples.keys())
    
    print('\n'.join(candidate_samples))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", type=str,
                        help="path to the folder containing samtools-generated depth files")
    parser.add_argument("--no-coverage", dest='cov', action='store_false', help="if set, do not filter by coverage")
    parser.set_defaults(cov = True)
    parser.add_argument("--fit-step", dest='fit_step', action='store_true', help="if set, filter by fitting a step function")
    parser.set_defaults(fit_step = False)
    parser.add_argument("-d", "--depth", type=int, help="depth coverage threshold", default=50)
    parser.add_argument("-b", "--breadth", type=float, help="breadth coverage threshold", default=0.8)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)

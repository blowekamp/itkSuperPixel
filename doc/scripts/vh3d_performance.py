#!/usr/bin/env python

from __future__ import print_function

descprtion = """
This script is for performance testing of the SLIC filter on the
Visible Human Male data set. Specifically, it run over a varied number
of threads to help determine scalability. The results a accumulated
into a pandas dataframe which is then written to a csv file
"""

import SimpleITK as sitk
from timeit import default_timer as timer
import timeit
import os
import pandas as pd
import argparse


parser =  argparse.ArgumentParser(description="Scalability analysis  for the SLIC image filter.")
parser.add_argument('inputImageFilename', help="An image file used for input")
parser.add_argument('-t','--threads', nargs='+', default=[88,66,44,32,24,16,8,4,2,1], help='<Required> Set flag', type=int)
parser.add_argument('-r','--repeat', type=int, default=1)

args = parser.parse_args()

#vhfilename_filefile='/scratch/blowekamp/vhm_lab.mha'
#inputFilename='/dev/shm/vhm_lab.mha'
#inputFilename='vhm_lab_100.mha'



baseFilename=os.path.basename(args.inputImageFilename)
baseFilename=os.path.splitext(baseFilename)[0]

print("Reading \"{0}\"...".format(args.inputImageFilename))
img = sitk.ReadImage(args.inputImageFilename)
print("DONE")


num_threads =  args.threads

timings1 = []

for np in num_threads:
    slic = sitk.SLICImageFilter()
    slic.SetSuperGridSize([50,50,50])
    slic.SetNumberOfThreads(np)
    slic.EnforceConnectivityOff()
    t = timeit.repeat(lambda: slic.Execute(img), repeat=args.repeat, number=1)

    print ("np: {0} time: {1}".format(np,t))
    timings1.append(min(t))


timings2 = []

for np in num_threads:
    slic = sitk.SLICImageFilter()
    slic.SetSuperGridSize([50,50,50])
    slic.SetNumberOfThreads(np)
    slic.EnforceConnectivityOn()
    start = timer()
    out = slic.Execute(img)

    t = timeit.repeat(lambda: slic.Execute(img), repeat=args.repeat, number=1)

    print ("np: {0} time: {1}".format(np,t))
    timings2.append(min(t))

df = pd.DataFrame.from_items( [ ("Number Of Threads", num_threads),
                                ("Without Connectivity", timings1),
                                ("With Connectivity", timings2) ] )
print(df)

df.to_csv("performance_{0}.csv".format(baseFilename),index=False)

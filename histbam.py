import matplotlib.pyplot as plt
import matplotlib
import argparse
matplotlib.use("Agg")
from datetime import datetime

IMGDIR = "/exports/sascstudent/samvank/images"

"""
Possible adjustment for logarithmic bins:
import numpy as np

def create_log_bins(data, num_bins):
    min_value = np.min(data)
    max_value = np.max(data)
    log_min = np.log10(min_value)
    log_max = np.log10(max_value)
    log_bins = np.logspace(log_min, log_max, num_bins)
    return log_bins

data = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
num_bins = 10

log_bins = create_log_bins(data, num_bins)
print(log_bins)
"""

def findMaxRange(filename):
    maxlen = 1
    fle = open(filename)
    while line:=fle.readline():
        if line[0:2] == "SN":
            if line.split("\t")[1] == "maximum length:":
                maxlen = int(line.split("\t")[2])
                break
    fle.close()
    return maxlen


def createBuckets(filename, bins, bucketLength):
    bucketSize = bins
    if bucketLength <= 0:
        bucketLength = findMaxRange(filename) / (bucketSize-1)
    buckets = [0] * bucketSize 
    ranges = []
    fle = open(filename)
    while line:=fle.readline():
        if line[0:2] == "RL":
            rlength = int(line.split("\t")[1])
            ramount = int(line.split("\t")[2])
            buckets[int((rlength) // bucketLength)] += ramount
    fle.close()
    for i in range(bucketSize):
        ranges.append(i*bucketLength)
    return buckets, ranges


def parseConfig():
    parser = argparse.ArgumentParser(prog="histvar",
                                    description="""
                                    Plots histogram of read-lengths in a samtools-stat file.
                                    ('samtools stat yourbamfile.bam > yourbamfile.stats') 
                                    """)
    parser.add_argument("-b", "--bins", action="store", nargs=1, type=int, default=[100],
                        help="amount of bins")
    parser.add_argument("-w", "--width", action="store", nargs=1, type=int, default=[0],
                        help="Supply bin-width (by default widths are calculated based on size-range)")
    parser.add_argument("-l", "--log", action="store_true",
                        help="Draw x-axis logarithmically.")
    parser.add_argument("file")
    args = parser.parse_args()
    return args


args = parseConfig()
buckets, ranges = createBuckets(args.file, args.bins[0], args.width[0])

plt.hist(ranges, args.bins[0],weights=buckets, edgecolor=None)
plt.xlabel("read Length (bp)")
plt.ylabel("amount of reads")
plt.yscale("log")
if (args.log):
    plt.xscale("log")

now = datetime.now()
filename =  args.file.split("/")[-1]
name =  IMGDIR + "/" + filename +  now.strftime("_%m-%d_%H:%M") + ".png"

plt.savefig(name, dpi=300)


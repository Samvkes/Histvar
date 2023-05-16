from cyvcf2 import VCF
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import numpy as np
import argparse
import gzip
from datetime import datetime

IMGDIR = "/exports/sascstudent/samvank/images"

def isIncorrect(event_range: list, only_short: bool, all_variants: bool) -> str:
    if len(set(event_range)) < 2:
        return "range values should be distinct"
    elif event_range[1] < event_range[0]:
        return "max_value of range should be larger than min_value"
    elif not only_short and not all_variants and event_range[0] > -51 and event_range[1] < 51:
        return "Given range will only show short events from a set that contains none, use -s or -a instead" 
    elif only_short and event_range[0] < -50 and event_range[1] > 50:
        print("Warning: given range will have no effect when only looking at short events")
        return "" 
    else:
        return ""


def readVarWithRange(varfile: str, event_range: list, only_short: bool, all_variants: bool) -> tuple:
    lenList = []
    rest = 0
    restLabel = "<50bp events"
    if all_variants:
        for variant in VCF(varfile):
            reflen = len(variant.REF)
            altlen = len(variant.ALT[0])
            if reflen - altlen > 0 and reflen - altlen <= event_range[1]:
                lenList.append(reflen - altlen)
            elif reflen - altlen < 0 and reflen - altlen >= event_range[0]:
                lenList.append(reflen - altlen)
    elif only_short:
        for variant in VCF(varfile):
            reflen = len(variant.REF)
            altlen = len(variant.ALT[0])
            if abs(reflen - altlen) < 50:
                if reflen - altlen > 0 and reflen - altlen <= event_range[1]:
                    lenList.append(reflen - altlen)
                elif reflen - altlen < 0 and reflen - altlen >= event_range[0]:
                    lenList.append(reflen - altlen)
            else:
                rest+= 1
        restLabel = ">49bp events"
    else:
        reallyLong = open("rlong.txt", "w")
        for variant in VCF(varfile):
            reflen = len(variant.REF)
            altlen = len(variant.ALT[0])
            if abs(reflen - altlen) > 49:
                le = reflen - altlen
                if le > 800000: 
                    reallyLong.write(f"ref: {variant.REF}, \n alt: {variant.ALT}")
                if reflen - altlen > 0 and reflen - altlen <= event_range[1]:
                    lenList.append(reflen - altlen)
                elif reflen - altlen < 0 and reflen - altlen >= event_range[0]:
                    lenList.append(reflen - altlen)
            else:
                rest+= 1
    return (lenList, rest, restLabel)


def readVar(varfile: str, only_short: bool, all_variants: bool) -> tuple:
    lenList = []
    rest = 0
    restLabel = "<50bp events"
    variants = VCF(varfile)
    if all_variants:
        for variant in VCF(varfile):
            reflen = len(variant.REF)
            altlen = len(variant.ALT[0])
            lenList.append(reflen - altlen)
    elif only_short:
        for variant in VCF(varfile):
            reflen = len(variant.REF)
            altlen = len(variant.ALT[0])
            if abs(reflen - altlen) < 50:
                lenList.append(reflen - altlen)
            else:
                rest+= 1
        restLabel = ">49bp events"
    else:
        for variant in VCF(varfile):
            reflen = len(variant.REF)
            altlen = len(variant.ALT[0])
            if abs(reflen - altlen) > 49:
                lenList.append(reflen - altlen)
            else:
                rest+= 1
    #else:
        #counter = 0
        #oldVariant = 0
        #while True:
            #counter += 1
            #try:
                #if variant:= next(variants):
                    #if counter % 10000 == 0: 
                        #print(variant.CHROM)
                        #if oldVariant == variant or variant.CHROM == "chr22":
                            #break
                        #oldVariant = variant
                         
                    #reflen = len(variant.REF)
                    #altlen = len(variant.ALT[0])
                    #if abs(reflen - altlen) > 49:
                        #lenList.append(reflen - altlen)
                    #else:
                        #rest+= 1
                #else:
                    #return (lenList, rest, restLabel)
            #except:
                #continue
    return (lenList, rest, restLabel)


def parseConfig():
    parser = argparse.ArgumentParser(prog="histvar",
                                    description="""
                                    Plots histogram of indel-length of long (>49bp) events in a given vcf. 
                                    A horizontal bar shows the amount of events not otherwise plotted
                                    (the amount of shorter (>49bp) events in the default graph, 
                                    the amount of longer (<50bp) events when -s is supplied, 
                                    and none when -a is supplied)
                                    """)
    parser.add_argument("file")
    parser.add_argument("-n", "--nomarker", action="store_false",
                        help="don't display unplotted-event-number marker")
    parser.add_argument("-i", "--image", action="store", nargs="?", metavar=("NAME"), const="noneSupplied", default="noneSupplied",
                        help="specify image path and name")
    parser.add_argument("-b", "--bins", action="store", nargs=1, type=int, default=[100],
                        help="amount of bins")
    parser.add_argument("-r", "--range", action="store", nargs=2, metavar=("MIN", "MAX"), type=int, default=[0,0],
                        help="set range of event-lengths to look at, min max (inclusive)")
    parser.add_argument("-l", "--linear", action="store_true",
                        help="show linear y-axis instead of logarithmic")
    parser.add_argument("-c", "--check", action="store_true",
                        help="check whether input is a valid vcf file")
    variant_types_group = parser.add_mutually_exclusive_group()
    variant_types_group.add_argument("-s", "--short", action="store_true",
                        help="plot histogram of short (<49bp) variants instead")
    variant_types_group.add_argument("-a", "--all", action="store_true",
                        help="plot histogram of all variants instead")
    args = parser.parse_args()
    return args


sns.set_theme()
args = parseConfig()


if args.check:
    counter = 0
    r = "\033[31m"
    g = "\033[32m"
    s = "\033[0m"
    try:
        for variant in VCF(args.file):
            counter += 1
            pass
        print(f"{g}Everything seems to be in order! No errors were found while checking variants.{s}")
    except Exception as e:
        print(f"{r}One or more variants seem to be corrupted / formatted wrong.")
        print(f"(You may want to run 'bcftools view {args.file} > /dev/null' for more details)")
        print(f"We're grabbing the offending line:{s}\n")
        lineCounter = 0
        headers = True
        colnames = ""
        lastLine = ""
        for line in gzip.open(args.file, "r"):
            if headers == True:
                decoded_line = line.decode("utf-8")
                if decoded_line != "##":
                    if decoded_line[0] == "#":
                        colnames = decoded_line.strip()
                    else:
                        headers = False
            else:
                lineCounter += 1
            if lineCounter == counter:
                print(f"\t\t\t{r}{colnames}\nLast correct line:{s}\t{lastLine}\n{r}First incorrect line:{s}\t{line.decode('utf-8')}")
                break
            else:
                lastLine = line.decode("utf-8")

    exit()


if args.range != [0,0]:
    if isIncorrect(args.range, args.short, args.all):
        print(isIncorrect(args.range, args.short, args.all))
        exit()
    else:
        lenList, rest, restLabel = readVarWithRange(args.file, args.range, args.short, args.all)
else:
    lenList, rest, restLabel = readVar(args.file, args.short, args.all)

if args.nomarker and not args.all:
    plt.axhline(rest, color=(0.8,0.2,0), alpha=0.7, label=f"amount of {restLabel}")
    plt.legend()

if not args.linear:
    plt.yscale("log")


ar = np.array(lenList)
plt.xlabel("indel Length (bp)")
plt.ylabel("amount of events")
sns.histplot(data=ar, bins=args.bins[0], )
if args.image != "noneSupplied":
    name = args.image
else:
    now = datetime.now()
    filename =  args.file.split("/")[-1]
    name =  IMGDIR + "/" + filename +  now.strftime("_%m-%d_%H:%M") + ".png"

plt.savefig(name, dpi=300)
print(f"Saved plot as {args.image}\n"
      f"Amount of plotted events: {len(lenList)} (largest in: {max(lenList)}bp, largest del: {min(lenList)}bp)\n"
      f"Amount of not-plotted events: {rest}")


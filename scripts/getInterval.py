import argparse
import sys
import textwrap
import os
parser = argparse.ArgumentParser(prog='fastaKit',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Description
    Get or Mask desired region of a fasta file
    For mask region in mask.bed: python3 getInterval.py -all -mask --cutfile mask.bed -fasta target.fasta > masked.fa
    For get a desired region: python3 getInterval.py -c -low ${begin} -up ${end} -fasta target.fasta > desired.fa
'''))
parser.add_argument('-cutted',  metavar='<file>', help='input cut file')
parser.add_argument('-fasta', metavar='<file>', help='input fasta file')
parser.add_argument('-o', metavar='<dir>', help='input out dir')
parser.add_argument('-c',  default=False, action='store_true', help='get Repeats')
parser.add_argument('-csv', metavar='<file>', help='input csv file')
parser.add_argument('-up', metavar='<INT>', type=int, help='upperbound')
parser.add_argument('-low', metavar='<INT>', type=int, help='lowerbound')
parser.add_argument('--cutfile', metavar='<FILE>', help='cut info')
parser.add_argument('-mask',  default=False, action='store_true', help='mask')
parser.add_argument('-all',  default=False, action='store_true', help='mask all region')

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

def readCut(fh):
    cutinfo = {}
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        lining = line.split("\t")
        name = lining[0]
        if name != "rdna" and name != "gap":
            cor = (int(lining[1]), int(lining[2]))
        else:
            cor = (0,0)
        if name in cutinfo:
            cutinfo[name].append(cor)
        else:
            cutinfo[name] = [cor]
    return cutinfo


def readCsv(fh):
    i = 0
    intervals = []
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        i = i + 1
        if i == 1:
            continue
        else:
            cutline = line.split(",")
            start = int(cutline[2])
            end = int(cutline[3])
            name = ">%s_%d-%d"%(cutline[1].replace("\"", ""), start, end)
            intervals.append([start,end,name])

    return intervals

def specialPre(seq,seqlenth):
    seqnew = []
    output = []
    n = len(seq)
    i = 0
    while i<n:
        if len(seq[i]) == seqlenth:
            seqnew.append(seq[i])
            i = i + 1
        elif len(seq[i]) > seqlenth:
            seqnew.append(seq[i][0:seqlenth])
            seq[i] = seq[i][seqlenth:]
        else:
            if i+1 >= n:
                seqnew.append(seq[i])
                break
            else:
                seq[i+1] = seq[i] + seq[i+1]
                i = i + 1
    output = seqnew
    return output


def readRe(fh):
    intervalcol = []
    interval = []
    start = 0
    end = 0
    count = 0
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n","")
        if line == "" and count == 0:
            #intervalcol.append([start, end])
            count = count + 1
            continue
        elif line != "" and not line.startswith("@"):
            lining = line.split("\t")
            start = int(lining[2])
            end = int(lining[3])
            count = 0
            if interval == []:
                interval = [start,end]
                continue
            elif interval[1]+ 200 >= start:
                interval[1] = end
                continue
            elif interval[1]+200 < start:
                intervalcol.append(interval)
                interval = [start,end]
                continue
    intervalcol.append(interval)
    return intervalcol

def readFasta(fh):
    seq = []
    info = ""
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if not seq == []:
                yield "".join(seq), info
            info = line.replace("\n","")
            seq.clear()
            continue
        else:
            seq.append(line.replace("\n",""))
    if info != "":
        yield "".join(seq), info

def cutSeq(interval, info, seq, check):
    for inter in interval:
        if inter == []:
            print("Empty File")
            break
        #if "=" in info:
         #   oldtitle = info.split("=")[0]
          #  iival = info.split("=")[1]
           # iistart = int(iival.split("-")[0]) + inter[0]
            #iiend = int(iival.split("-")[0]) + inter[1]
            #title = "%s=%d-%d" % (oldtitle, iistart, iiend)
        #else:
            #title = "%s=%d-%d" % (info, inter[0], inter[1])
        len = inter[1] - inter[0]
        #if not args.c and not check:
            #if len <3000:
            #    continue
        #if check:
            #title = inter[2]
        title = info.split(" ")[0]
        title = "%s:%d-%d" % (title, inter[0], inter[1])
        newseq = seq[inter[0]:inter[1]]
        #newnewseq = cutseq(newseq)
        yield title, newseq

def cutseq(seq):
    newseq = []
    reminder = len(seq) % 100
    iter = len(seq) // 100
    i = 0
    while i<iter:
        newseq.append(seq[i*100:(i+1)*100])
        i = i + 1
    newseq.append(seq[iter*100:])
    return "\n".join(newseq)

def mask(interval, info, seq):
    for inter in interval:
        if inter == []:
            print("Empty")
            break
        len = inter[1] - inter[0]
        title = info.split(" ")[0]
        masking = "N" * len
        newseq = seq[0:inter[0]] + masking + seq[inter[1]:]
        title = title + "_mask:%d-%d"%(inter[0],inter[1])
        yield title, newseq

def maskall(interval, info, seq):
    newseq = ""
    end = 0
    #print(interval)
    interval.sort()
    #print(interval)
    i = 0
    while i < len(interval):
        inter = interval[i]
        if inter == []:
            print("Empty")
            break
        length = inter[1] - inter[0]
        #print(inter)
        masking = "N" * length
        #print(length)
        #print(newseq)
        newseq = newseq + seq[end:inter[0]] + masking
        end = inter[1]
        i = i + 1
    if inter[1] < len(seq):
        newseq = newseq + seq[end:]
    title = info.split(" ")[0]
    title = title + " mask_via_%s" % (args.cutfile)
    return title, newseq


if __name__ == "__main__":
    if not (args.c or args.mask):
        if not os.path.exists(args.o):
            os.mkdir(args.o)
        dirpath = args.o
        if args.cutted:
            with open(args.cutted, "r") as fh:
                interval = readRe(fh)
    with open(args.fasta, "r") as f:
        if args.cutfile:
            with open(args.cutfile,"r") as fh:
                cutinfo = readCut(fh)
        #print(cutinfo)
        for seq,info in readFasta(f):
            if args.c:
                if args.cutfile:
                    seqname = info.replace(">","").split(" ")[0]
                    if seqname in cutinfo:
                        interval = cutinfo[seqname]
                        for title, newseq in cutSeq(interval, info, seq, False):
                            print(title)
                            NEWSEQ = newseq.upper()
                            print(NEWSEQ)
                else:
                    interval = [[args.low,args.up]]
                    for title, newseq in cutSeq(interval, info, seq, False):
                        print(title)
                        NEWSEQ = newseq.upper()
                        print(NEWSEQ)
            elif args.cutted:
                for title, newseq in cutSeq(interval,info,seq, False):
                    pathname = "%s/%s.fasta" % (dirpath, title.split(">")[1])
                    with open(pathname,"a") as out:
                        out.write(title)
                        out.write("\n")
                        out.write(cutseq(newseq))
            elif args.csv:
                with open(args.csv, "r") as fh:
                    interval = readCsv(fh)
                    pathname = "%s/HOR_%s.fasta" % (args.o, info.replace(">",""))
                    for title, newseq in cutSeq(interval,info, seq, True):
                        with open(pathname, "a") as out:
                            out.write(title)
                            out.write("\n")
                            out.write(newseq)
                            out.write("\n")
            elif args.mask:
                #print("1")
                if not args.cutfile:
                    length = args.up-args.low
                    masking = "N" * length
                    sequence = seq[0:args.low] + masking + seq[args.up:]
                    name = info + "_mask:%d-%d"%(args.low,args.up)
                    print(name)
                    print(cutseq(sequence))
                else:
                    seqname = info.replace(">","").split(" ")[0]
                    #print("2")
                    #print(seqname)
                    #print(cutinfo)
                    if seqname in cutinfo:
                        #print("3")
                        interval = cutinfo[seqname]
                        #print(interval)
                        if args.all:
                            title, newseq = maskall(interval, info, seq)
                            print(title)
                            print(newseq.upper())
                        else:
                            for title, newseq in mask(interval, info, seq):
                                print(title)
                                print(newseq)
                    else:
                        print(info)
                        print(seq)





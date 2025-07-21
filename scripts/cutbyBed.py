import argparse
import sys
import textwrap
import os
import gzip
from Bio.Seq import Seq

parser = argparse.ArgumentParser(prog='fastaKit',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Description
    get subregion of a fasta file from a given bed file.
    usage python3 cutbyBed.py -fasta target.fa -bed region.bed > region.fa
'''))

parser.add_argument('-fasta', metavar='<file>', help='input fasta file')
parser.add_argument('-bed', metavar='<file>', help='input bed file')
parser.add_argument('-addinfo', default=False, action='store_true', help='add 5th column in the seuqence name')
parser.add_argument('-infoonly', default=False, action='store_true', help='set sequence name as 5th column')
parser.add_argument('-dic', default=False, action='store_true', help='make reverse complement when the 6th column of the bed file is -')

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

def readBed(fh):
    cutinfo = {}
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        lining = line.split("\t")
        name = lining[0]
        if len(lining) >= 4:
            type = lining[3]
        else:
            type = "na"
        if len(lining) <= 5:
            dic = "na"
        else: 
            dic = lining[5]
        cor = (int(lining[1]), int(lining[2]),dic,type)
        if name in cutinfo:
            cutinfo[name].append(cor)
        else:
            cutinfo[name] = [cor]

    return cutinfo

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
        yield "".join(seq), info.split(" ")[0]

def cutSeq(interval, info, seq, check):
    for inter in interval:
        if inter == []:
            print("Empty File")
            break
        len = inter[1] - inter[0]
        #if not args.c and not check:
            #if len <3000:
            #    continue
        #if check:
            #title = inter[2]
        title = info.split(" ")[0]
        #print(inter)
        if args.addinfo:
            title = "%s:%d-%d_{%s}" % (title, inter[0], inter[1],inter[3])
        elif args.infoonly:
            title = ">%s %s:%d-%d" % (inter[3],title.replace(">",""), inter[0], inter[1])
        else:
            title = "%s:%d-%d" % (title, inter[0], inter[1])
        newseq = seq[inter[0]:inter[1]]
        #print(inter[2])
        if args.dic and inter[2] == "-":
            #print(inter[2])
            seqre = str(Seq(newseq).reverse_complement())
            title =title + "_R"
            yield title, seqre
        else:
        #newnewseq = cutseq(newseq)
            yield title, newseq

def main():
    with open(args.bed, "r") as fh:
        cutinfo = readBed(fh)
    with open(args.fasta, "r") as f:
        for seq,info in readFasta(f):
            seqname = info.replace(">","").split(" ")[0]
            if seqname in cutinfo:
                interval = cutinfo[seqname]
                for title, newseq in cutSeq(interval, info, seq, False):
                    print(title)
                    NEWSEQ = newseq.upper()
                    print(NEWSEQ)

if __name__ == "__main__":
    main()

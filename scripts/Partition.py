import argparse
import sys
import textwrap
import os

parser = argparse.ArgumentParser(prog='boxplt',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Description
    get duplication id
'''))

parser.add_argument('-rel', metavar='<file>', help='input file')


if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

def checker(idn,repcol):
    for ids in repcol.keys():
        partition = repcol[ids]
        if idn in partition:
            #repcol[ids].append(idn)
            return ids, True
    return id,False

with open(args.rel, 'r') as fh:
    repcol = {}
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        lining = line.split("\t")
        id1 = lining[0]
        id2 = lining[1]
        id1_ids,id1_check = checker(id1,repcol)
        id2_ids,id2_check = checker(id2,repcol)
        #print(repcol)
        if id1_check and id2_check:
            #print("1")
            continue
        elif id1_check:
            repcol[id1_ids].append(id2)
            #print("2")
            continue
        elif id2_check:
            repcol[id2_ids].append(id1)
            #print("3")
            continue
        else:
            if id1 == id2:
                repcol[id1] = [id1]
            else:
                repcol[id1] = [id1 ,id2]
            #print("4")

for ids in repcol:
    partition = repcol[ids]
    print("\t".join(partition))
    

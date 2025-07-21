import sys,os
import re
f1 = open(sys.argv[1],'r')
dicts = {}
for line in f1:
    chr1,start1,end1,dirs,geneid,chrs,start,end = line.strip().split()
    key = "%s:%s-%s"%(chrs,start,end)
    if key not in dicts:
        dicts[key] = {}
    if geneid not in dicts[key]:
        dicts[key][geneid] = [0,[]]
    dicts[key][geneid][0] += 1
    dicts[key][geneid][1].append(dirs)
f1.close()

for key in dicts:
    localdict = dicts[key]
    if "Foxk2" in localdict and "OR4" in localdict:
        Model = "M1"
    elif "OR4" in localdict and "Foxk2" not in localdict:
        Model = "M2"
    elif "OR4" not in localdict and "Foxk2" in localdict:
        Model = "M3"
    else:
        Model = "Others"
    srtk = ""
    for k,v in localdict.items():
        srtk += "%s:%s:%s; "%(k,v[0],','.join(v[1]))
    print('%s\t%s\t%s'%(key,Model,srtk))

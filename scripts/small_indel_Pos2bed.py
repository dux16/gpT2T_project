import sys

if len(sys.argv) < 2:
    sys.exit("python3 {} *.indel".format(sys.argv[0]))
def parse_pos(string):
    if "-" in string:
        # length > 1
        pos1 = int(string.split("-")[0]) - 1   # bed format
        pos2 = int(string.split("-")[1])
    else:
        # length = 1
        pos2 = int(string) # bed format, 1-POSITION 
        pos1 = pos2 - 1 # bed format, 0-POSITION 


    return (pos1, pos2)

with open(sys.argv[1], 'r') as fh:
    for i in fh:
        if i.startswith("#"):
            continue
        tmp = i.strip().split("\t")
        if tmp[4] == 'ins':
            ref_start = int(tmp[2].replace("_", "")) - 1 # bed format
            ref_end = ref_start
            (que_start, que_end) = parse_pos(tmp[3])
            print(f"{tmp[0]}\t{ref_start}\t{ref_end}\t{tmp[1]}\t{que_start}\t{que_end}\tins")
        else:
            que_start = int(tmp[3].replace("_", "")) - 1 # bed format
            que_end = que_start
            (ref_start, ref_end) = parse_pos(tmp[2])
            print(f"{tmp[0]}\t{ref_start}\t{ref_end}\t{tmp[1]}\t{que_start}\t{que_end}\tdel")

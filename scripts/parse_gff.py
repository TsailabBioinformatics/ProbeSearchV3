import sys

clean = open("./clean.table", "wt")

def parse_gff(gff):
    f = open(gff)
    for line in f:
        line = line.split("\t")
        line[3] = line[3][(line[3].find("=") + 1):find_nth(line[3], ".", 2)]
        # line[3] = line[3][:16]
        clean.write("\t".join(line))
        clean.write("\n")
        

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start


gff = sys.argv[1]
parse_gff(gff)
import sys

clean = open("./clean.table", "wt")

def parse_gff(gff):
    count = 0
    f = open(gff)
    for line in f:
        line = line.split("\t")
        if (len(line) == 4):
            count = count + 1
            #line[0] = line[0].replace("Chr", "A")
            line[3] = line[3][(line[3].find("=") + 1):find_nth(line[3], ".", 2)]
            #line[3] = line[3][:line[3].find(";")]
            clean.write("\t".join(line))
            clean.write("\n")
    print(count)
            

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start


gff = sys.argv[1]
parse_gff(gff)
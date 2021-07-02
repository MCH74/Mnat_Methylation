#! /usr/bin/python3

import sys
counter = 0
with open(sys.argv[1],"r") as f:
    for line in f:
        if ">" not in line:
            line = line.strip()
            counter += line.upper().count("CG")
print(counter)


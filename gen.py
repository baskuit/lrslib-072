"""
Prints random mxn constant-sum (=den) matrix. Lrsnash format.
"""

import sys
import random

if len(sys.argv) != 4:
    print("Usage: python script.py rows cols den")
    sys.exit(1)

# Extract command-line arguments
rows = int(sys.argv[1])
cols = int(sys.argv[2])
den = int(sys.argv[3])

print(rows, cols, '\n')

array = [[0 for __ in range(cols)] for _ in range(rows)]

for r in range(rows):
    for c in range(cols):
        array[r][c] = random.randint(0, den - 1)

for p in range(2):
    for r in range(rows):
        s = ""
        for c in range(cols):
            x = array[r][c]
            s += str(x) if p == 0 else str(den - x)
            s += " "
        print(s)
    print()


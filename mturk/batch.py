from __future__ import print_function, division, unicode_literals
import argparse
import csv
import json
import random
import os
import sys
import traceback

from typing import List

def write_batches(batches: List, filename: str):
    with open(filename, 'w') as file:
        for b in batches:
            file.write(json.dumps(b) + "\n")

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def main(options):
    input = options.input
    if input.name.endswith('.csv'):
        entries = list(csv.DictReader(input))
    elif input.name.endswith('tab') or input.name.endswith('tsv'):
        entries = list(csv.DictReader(input, delimiter='\t'))
    else:
        # Assume one line per entry
        entries = [line.strip() for line in input.readlines() if len(line.strip()) > 0]

    if options.shuffle:
        random.shuffle(entries)
    batches = chunks(entries, options.n)
    write_batches(batches, options.output)


def parse_arguments():
    parser = argparse.ArgumentParser(usage="Batch items and outputs them")
    parser.add_argument('-n', type=int, help="Batch size")
    parser.add_argument('--shuffle', action="store_true", help="Shuffle input")
    parser.add_argument('-i', '--input', type=argparse.FileType("r"), help="Input filename")
    parser.add_argument('-o', '--output', help="Output filename")
    options = parser.parse_args()
    return options

if __name__ == "__main__":
    # Parse options
    options = parse_arguments()

    # Run the program
    code = 0
    try:
        main(options)
    except:
        traceback.print_exc()
        code = -1
    sys.exit(code)

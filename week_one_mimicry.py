#!./venv/bin/python

from bio_functions import *
import sys

print(sys.argv[0])

data = dir_h5_to_np("data")

for d in data:
    print(f"{d['label']} {type(d['keys'])} {type(d['data'])}")

#!./venv/bin/python

from bio_functions import *
import sys

print(sys.argv)

if not os.path.isdir(sys.argv[1]):
    print(f"Directory |{sys.argv[1]}| not found")
    sys.exit()

data = dir_h5_to_np(sys.argv[1])

for d in data:
    print(f"{d['label']} {type(d['keys'])} {type(d['data'])}")

#!/usr/bin/env python

import matplotlib
matplotlib.rc("figure", max_open_warning=0)
from matplotlib import pyplot as plt

import capnp
import hit_capnp

def main():
    f = open("data/testout.hits")
    hits = hit_capnp.Hit.read_multiple(f)
    for hit in hits:
        print(hit)
    
if __name__ == "__main__":
    main()

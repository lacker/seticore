#!/usr/bin/env python

import ctypes
import os
import sys

assert __name__ == "__main__"

DIR = os.path.dirname(os.path.realpath(__file__))
handle = ctypes.CDLL(f"{DIR}/hello.so")
handle.hello_world.argtypes = [ctypes.c_int]

handle.hello_world(1337)

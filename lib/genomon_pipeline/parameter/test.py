#! /usr/bin/env python

from genomon_param import genomon_param
import test_func

genomon_param.read("param.cfg")

print genomon_param.get("bwa", "option")

test_func.test_func()




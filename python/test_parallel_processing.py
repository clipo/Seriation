__author__ = 'clipo'

import multiprocessing

try:
    cpus = multiprocessing.cpu_count()
except NotImplementedError:
    cpus = 2   # arbitrary default

print "number of cpus" , cpus

def square(n):
    return n * n

pool = multiprocessing.Pool(processes=cpus)
print pool.map(square, xrange(1000))
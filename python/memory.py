__author__ = 'carllipo'

import argparse
import sys
import resource

class Memory():
    def __init__(self):
        self.previousMemory=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
        return

    def memory(self):
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

    def memoryChange(self):
        return (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000) - self.previousMemory


def main():
    mem=Memory()
    dict = {}
    for i in range(0,1000):
        dict[i] = 'abcdefg'
        print mem.memory()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='track memory')
    parser.add_argument('--test')
    try:
        args = vars(parser.parse_args())
    except IOError, msg:
        parser.error(str(msg))
        sys.exit()
    if args['test'] is not None:
        main()
'''
This is the main file that will be used to run the black hole operations.

Neetre 2024
'''

import argparse

from ev_bh import event_horizon
from bh_lensing import lensing_bh
from bh_mergin import merge_bh
from binary_stars import binary_stars


def args_parser():
    parser = argparse.ArgumentParser(description='Black Hole Operations')
    parser.add_argument('-eh', '--event_horizon', action='store_true', help='Calculate the event horizon of a black hole')
    parser.add_argument('-lbh', '--lensing_bh', action='store_true', help='Calculate the lensing of a black hole')
    parser.add_argument('-mbh', '--merge_bh', action='store_true', help='Calculate the merging of two black holes')
    parser.add_argument('-bs', '--binary_stars', action='store_true', help='Calculate the motion of binary stars')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    return parser.parse_args()


def main():
    args = args_parser()
    print(args)
    
    if args.event_horizon:
        event_horizon(False, False)
    
    if args.lensing_bh:
        lensing_bh()

    if args.merge_bh:
        merge_bh()
    
    if args.binary_stars:
        binary_stars()
    
    
if __name__ == "__main__":
    main()
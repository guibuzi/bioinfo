import os, sys, argparse

def cmd():
    args = argparse.ArgumentParser(description='Recombinntion analysis.', epilog='Jinfeng Zeng')
    args.add_argument('--version', action='version', version='version 1.0')
    aa = args.parse_args()
    return aa
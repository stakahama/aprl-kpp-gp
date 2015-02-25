import argparse

parser = argparse.ArgumentParser(description='build kpp for gas-phase simulation')
parser.add_argument('ROOT',type=str)
parser.add_argument('KPPC',type=str)
parser.add_argument('CPATH',type=str)
parser.add_argument('--skipbuild', dest='skipbuild', action='store_true')
parser.set_defaults(skipbuild=False)
args = dict(vars(parser.parse_args()))

print args


import os
import sys         # FB: for command line arguments
dd = lambda x: os.path.dirname(os.path.dirname(x))
sys.path.append(os.path.join(dd(__file__),'lib'))

## user-defined
from matmul import matmul
from simpolmod import Simpolclass





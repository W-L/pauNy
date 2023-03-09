import sys
from pathlib import Path
prog_path = Path(sys.argv[0]).parent
sys.path.insert(0, f"{prog_path}/pauNy")
sys.path.insert(0, "./pauNy")
from lib import *
from plot import *

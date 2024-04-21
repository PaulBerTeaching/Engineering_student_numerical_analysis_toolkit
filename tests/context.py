import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append("../root_finding_problem")

import root_finding_problem.root_finding_methods as tf
from numpy import polynomial as poly
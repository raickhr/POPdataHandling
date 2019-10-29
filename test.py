import glob
import sys
import os
from glob import glob
import pandas as pd
import numpy as np

from readWrite import *
from gridModule import *
from operators import *
from getProdFields import writeWithProdFields

POP_outFileList = glob(POP_outpath+'*')
for fileName in POP_outFileList:
    writeWithProdFields(fileName,PRD_outpath)

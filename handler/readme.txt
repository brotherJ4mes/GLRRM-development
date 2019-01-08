# Downloaded by ZAM on 12/4/18 from https://github.com/tim-hunter/GLRRM 
#
# databank.py
# databank_io.py
# databank_util.py

IMPORTANT: From the three files downloaded from online, ZAM made three changes from

import databank as databank
import databank_util as util

to 

import handler.databank as databank
import handler.databank_util as util

BECAUSE: ZAM placed the files into a directory called "handler" and is calling them from a script outside of this folder.

# ZM NOTE 11/28 databank_util.py had original line as databank as bank




#/bin/python
# UNOFFICIAL SCRIPT - just for sharing between Yin,Tim,Zoe and Lauren so
# that we are all on the same page. One year later it is expected that this 
# script will drastically change and clean up...

#
#----------------------------------------------------------------
#  This is a really basic example of how the overall GLRRM might
#  be structured to use the databank repository.  The actual
#  implementation of GLRRM will be very different, once it is
#  developed.  This is mainly intended as a demonstration of how
#  to use the repository, while also illustrating the basic structure
#  that Tim Hunter had in mind when he developed the design documents.
#-----------------------------------------------------------------
#  This wrapper script is just a work in progress for testing functionalities of
#  many things - data handler, modules, etc.

# ----------------- import modules from python -------------------
import sys
sys.path
# ... zoe isn't even clear what the above does ... 

# ------------ import created module/package: handler ----------------
# three files for the databank package in folder and readme file


from handler import databank as db
# from handler import databank_util as db_util
from handler import databank_io as db_io


# ------ import created module/package: superior & middle & ontario -----------

# sys.path.insert(0,'/superior')
# sys.path.insert(0,'/middle')
# sys.path.insert(0,'/ontario)

from superior import routesuperior as sp
#import superior.superior_util as sp_util # tools and functions unique to superior
from middle import routemiddle as ml
from middle import routemiddle_util as ml_util # tools and functions unique to middle
# import ontario.routeontario as on
# import ontario.routeontario_util as on_util
 
help(ml.routemiddle)
help(ml.SolveIterative)
help(sp.routesuperior)
#------------------------------------------------------------------------------
# TODO: Create Configuration Files
#------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Create Vaults for Programs (Data and Record)
#-------------------------------------------------------------------------------
# !!! NOTE: if you change this below - you will need to update the config file to 
# match ... 

# ONLY TRY DAILY AND MONTHLY AT THE MOMENT

# kstep = 'dy'
# kstep = 'wk'
# kstep = 'qm'
kstep = 'mn'

folderin = 'data-in/' + kstep

from prep_vaults import create_vault_record, create_vault_all

# See "create_vault_all" for details about the .txt file stuctues expected
VaultRecord = create_vault_record(folderin)
VaultData = create_vault_all(folderin)

# Review what is in the Vaults
VaultRecord.printVault()
VaultData.printVault()


#---------------------------------------------------------------------------------
# ...temporary "hack" of pushing sup & ont level and flows into VaultData ...
#--------------------------------------------------------------------------------
# ...this step of pushing the actual sup and ont level and flow is so that
# in the dummy regulation plans they can be pulled down and pushed back in ...

# withdraw the superior and ontario outflow and levels (would be module output)
sup = VaultRecord.withdraw(kind='mlv', units='m', intvl=kstep, loc='sup')
VaultData.deposit(sup)
sup = VaultRecord.withdraw(kind='flw', units='cms', intvl=kstep, loc='smr')
VaultData.deposit(sup)

ont = VaultRecord.withdraw(kind='mlv', units='m', intvl=kstep, loc='ont')
VaultData.deposit(ont)
ont = VaultRecord.withdraw(kind='flw', units='cms', intvl=kstep, loc='stl')
VaultData.deposit(ont)

# review some information for data series ... 
sup.getOneLineSummary()
ont.dataVals[0:4]

# testing the interval retrieval
# testinvl = VaultData.withdraw(kind='mlv', units='m', intvl='mn', loc='sup').dataVals
# z = VaultRecord.withdraw(kind='flw', units='cms', intvl='dy', loc='stl')

#---------------------------------------------------------------------------------
# Read Module Configuration Files
#--------------------------------------------------------------------------------

# sp_file = 'config/config_superior.cfg'
ml_file = 'config/config_middle.cfg'
# on_file = 'config/config_ontario.cfg'

# sp_info = sp_util.sp_config_info(sp_file)
ml_info = ml_util.ml_config_info(ml_file)
# on_info = on_util.on_config_info(on_file)

# create a function in ml_info to print out the entire contents...to view

#------------------------------------------------------------------------------
# Simulation Start and End Date: Recall data read in could be more than needed
#------------------------------------------------------------------------------

# of format date (should actually be read from the master config...)
sim_start = ml_info.sdate
sim_end   = ml_info.edate

print(sim_start)
print(sim_end)

# shall be of format
# yyyy-mm-dd
# yyyy-mm-dd


#------------------------------------------------------------------------------
# Run The Models: Superior, Middle, Ontario
#------------------------------------------------------------------------------

[sp_run_return, spdf_header, spdf_file, storsp] = sp.routesuperior(VaultData, sim_start, sim_end, ml_info)
[ml_run_return, mldf_header, mldf_file, storml] = ml.routemiddle(VaultData, sim_start, sim_end, ml_info)
# [on_run_return, ondf_header, ondf_file, storon] = on.routeontario(VaultData, sim_start, sim_end, on_info)


# what is returned? 
# ...4 lists ... see code for more detail

#---------------------------------------------------------------------------------
# POST PROCESS: WRITE TEXT FILE OUTPUT
#--------------------------------------------------------------------------------

print(kstep)

# Withdraw level dataseries kind for outputing text files...mhu, stc, eri
sim_res_mlv_mhu = VaultData.withdraw(kind='mlv', units='m', intvl=kstep, loc='mhu', first=sim_start, last=sim_end)
sim_res_mlv_stc = VaultData.withdraw(kind='mlv', units='m', intvl=kstep, loc='stc', first=sim_start, last=sim_end)
sim_res_mlv_eri = VaultData.withdraw(kind='mlv', units='m', intvl=kstep, loc='eri', first=sim_start, last=sim_end)
	
# Withdraw flow dataseries kind for outputing text files...st.clair, detroit, niagara
sim_res_flw_scr = VaultData.withdraw(kind='flow', units='cms', intvl=kstep, loc='scr', first=sim_start, last=sim_end)
sim_res_flw_det = VaultData.withdraw(kind='flow', units='cms', intvl=kstep, loc='det', first=sim_start, last=sim_end)
sim_res_flw_nia = VaultData.withdraw(kind='flow', units='cms', intvl=kstep, loc='nia', first=sim_start, last=sim_end)

# write the results to files in the data-out folder
name_fold = 'data-out' + '/' + kstep + '/'

# write the files...
db_io.write_file(name_fold + 'middle_mlv_mhu.txt',"column",sim_res_mlv_mhu,overwrite=True)
db_io.write_file(name_fold + 'middle_mlv_stc.txt',"column",sim_res_mlv_stc,overwrite=True)
db_io.write_file(name_fold + 'middle_mlv_eri.txt',"column",sim_res_mlv_eri,overwrite=True)

db_io.write_file(name_fold + 'middle_flw_mhu.txt',"column",sim_res_flw_scr,overwrite=True)
db_io.write_file(name_fold + 'middle_flw_stc.txt',"column",sim_res_flw_det,overwrite=True)
db_io.write_file(name_fold + 'middle_flw_eri.txt',"column",sim_res_flw_nia,overwrite=True)

#----------------------------------------------------------------------------
# To Do - Within Middle Lakes - Write Output and Computation/Checks
#----------------------------------------------------------------------------
# this here is scrappy right now --- just first stab at writing output
# https://realpython.com/python-csv/#writing-csv-file-from-a-dictionary-with-csv
import numpy as np
import pandas as pd

# turn the data into a dataframe
mldf = pd.DataFrame(np.transpose(ml_run_return))

# write the file
mdlf_path = name_fold + mldf_file
mldf.to_csv(mdlf_path,index=False,header=mldf_header) 




#------------------------------------------------------------------------------
# POST PROCESS: PLOTS 
#------------------------------------------------------------------------------
from post_plots import make_hydrologic_plots
from post_plots import make_compare_plots

# Most of the data is in the Vault
# Storage Values are currently in the ml_run_return

make_hydrologic_plots(VaultData,storml,storsp,ml_info)
make_compare_plots(VaultData,VaultRecord,ml_info)


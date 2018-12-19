# SUPER UNNOFICIAL SCRIPT - to be edited and updated later

# ----------------------------------
# Function to read in the the data files using the data handler
# ----------------------------------

import handler.databank as db
# import handler.databank_util as db_util
import handler.databank_io as db_io

from os import listdir

#---------------------------------------------------------------------------------
# Create Vault for time series of observed time series ... 
#--------------------------------------------------------------------------------
def create_vault_record(folderdat):
    """
    INPUT:
        folder directory of input files
    
    ABOUT:
        will only read in the file if prefix is "record"
        
    """
    # create empty Vault
    VaultRecord = db.DataVault()
    
    # create list of files in the directory
    filelist = listdir(folderdat)
    # loop through files and only read in those with prefix "record"
    for n in range(0,len(filelist)):
        
        namefile = folderdat + '/' + filelist[n]
        print(namefile)

        if filelist[n][0:6]=='record':
           print(filelist[n][0:6])
           record_dat = db_io.read_file(filename=namefile)
           VaultRecord.deposit(record_dat)
     
    # output shall be VaultRecord
    return VaultRecord


#---------------------------------------------------------------------------------
# Create - Vault for  DAILY, WEEKLY, QUARTER-MONTHLY, MONTHLY
#--------------------------------------------------------------------------------

def create_vault_all(folderdat):
    """
    INPUT:
        folder directory of input files
    
    ABOUT:
        will only read in the file if prefix is NOT "record"
        
    """
    
    # create empty Vault
    Vault = db.DataVault()
    
    # create list of files in the directory
    filelist = listdir(folderdat)
    
    # loop through files and only read in those without prefix "record"
    for n in range(0,len(filelist)):
        namefile = folderdat + '/' + filelist[n]
        print(namefile)
        if filelist[n][0:6]!='record':
           dat = db_io.read_file(filename=namefile)
           Vault.deposit(dat)
     
    # output shall be Vault
    return Vault


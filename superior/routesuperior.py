#/bin/python
# from handler.databank import "?" ????
import handler.databank as db
#----------------------------------------------------------------

#-----------------------------------------------------------------

def routesuperior(VaultData, SimulationStartDate, SimulationEndDate, SuperiorInfo):
    """
    Lake Superior Routing: Solve for period lake level and flow from the start
    to end date
    
    ABOUT: 
    Function solves the "hydrologic equations" for all periods (routing steps)
    by breaking the problem into incremental steps and using regulation logic. 
    
    INPUTS: 4 Total (all required)
    1 object of DataVault type (see handler/databank for more detail)
            - contains time series to be withdrawn
    2 variables of date type:  SimulationStartDate, SimulationEndDate
    1 object of sp_config_info: (see superior_util for more detail)
            - contains crtieria for superior routing and settings
    
    OUTPUT: Results are pushed to vault during function/program.
           at the moment nothing is actually returned...
    """
    #-----------------------------------------------------------------------
    # CHECK INPUTS WERE PROVIDED [Basic Check - ATM no check for data type and dates]
    #-----------------------------------------------------------------------
    if not VaultData:
        raise Exception('Lake Superior Program missing VaultData')
    if not SimulationStartDate:
        raise Exception('Lake Superior Program missing SimulationStartDate')
    if not SimulationEndDate:
        raise Exception('Lake Superior Program missing SimulationEndDate')
    if not SuperiorInfo:
        raise Exception('Lake Superior Program missing SuperiorInfo')
	
	
    #-----------------------------------------------------------------------
    # UNPACK VARIABLES NEEDED FOR SOLUTION CALCULATIONS
    #-----------------------------------------------------------------------
    
    # "Unpack" the 2 dates
    startdate = SimulationStartDate
    finaldate = SimulationEndDate
    
    # define the resolution
    #key_intvl = SuperiorInfo.keyintvl
    
    key_intvl = SuperiorInfo.keyintvl
    oneday = SuperiorInfo.timedelta_one_day
    
    # define the constants (float)
    zero = 0.0
    half = 0.5
    secinday = 86400
    
    solvetestcms = 75
    start_lev_sup = 183.58 # change later to read off the config file
    spverbosity = 6
    spincnumber = SuperiorInfo.mlincrements    # change later to spincrements 
    spsolutionid = 1 # 1 = Force , 2 = ???
    
    Asp = 82100000000.00 # m2
    
    try:
        # Withdraw Net Basin Supplies   
        list_nbs_sup = VaultData.withdraw(kind='nbs', units='cms', 
            intvl=key_intvl, loc='sup', first=startdate, last=finaldate).dataVals
        list_len = len(list_nbs_sup)
    except:
        raise Exception('Lake Superior: Missing Net Basin Supplies')
        
    try:
        # Withdraw Diversion   
        list_div_sup = VaultData.withdraw(kind='flw', units='cms', 
            intvl=key_intvl, loc='oll', first=startdate, last=finaldate).dataVals
    except:
        list_div_sup = [150.00]*list_len
        raise Exception('Lake Superior: Missing Long-Lac Ogoki Diversion')
       
    try:
        # Withdraw Ice and Weed Retardation
        list_icw_sup = VaultData.withdraw(kind='icw', units='cms',
            intvl=key_intvl, loc='smr', first=startdate, last=finaldate).dataVals
    except:
        raise Exception('Lake Superior: Missing Ice and Weed Retardation Values')
    
    try:
        # Withdraw Consumptive Use
        list_csu_sup = VaultData.withdraw(kind='con', units='cms',
            intvl=key_intvl, loc='sup', first=startdate, last=finaldate).dataVals       
    
    except:
        list_csu_sup = [0]*list_len
        print('Lake Superior: Consumptive Use Missing from Vault, set to zeros')
    
    try:
        # Withdraw Groundwater
        list_gdw_sup = VaultData.withdraw(kind='flw', units='cms',
            intvl=key_intvl, loc='sup', first=startdate, last=finaldate).dataVals       
    
    except:
        list_gdw_sup = [0]*list_len
        print('Lake Superior: Consumptive Use Missing from Vault, set to zeros')    

    # COMPLETED: UNPACKING VARIABLES NEEDED FOR SOLUTION CALCULATIONS
    
    #----------------------------
    # dummy flows - reading in and writing out
    #----------------------------
    try:
        # Withdraw the forced level
        list_mlv_sup = VaultData.withdraw(kind='mlv', units='m',
            intvl=key_intvl, loc='sup', first=startdate, last=finaldate).dataVals
        list_flo_smr = VaultData.withdraw(kind='flw', units='cms',
            intvl=key_intvl, loc='smr', first=startdate, last=finaldate).dataVals   
    
    except:
        list_mlv_sup = [None]*list_len
        list_flo_smr = [None]*list_len
        print('Lake Superior: No Forced Flows Provided')  
    
    print(list_mlv_sup[0:5])
    print(list_flo_smr[0:5])
    #-----------------------------------------------------------------------
    # BEGIN LAKE SUPERIOR ROUTING PROGRAM
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    # INITIALIZATION
    #-----------------------------------------------------------------------
    print('The Start Date: ', startdate)
    print('The End Date: ', finaldate)

    # Create empty lists    
    #list_mlv_sup = [None]*list_len
    #list_flo_sup = [None]*list_len

    list_str_sup = [None]*list_len   
    
    list_dt_bop = [None]*list_len
    list_dt_eop = [None]*list_len
    list_dt_day = [None]*list_len
    list_dt_sec = [None]*list_len
    
    # Initialize BOP levels as starting levels for first period
    spbop = start_lev_sup

    # Initiate the routing date
    routedateBOP = startdate
    
    # Number of routing periods to go through
    routing_steps = list_len
    print('There are ', routing_steps,' routing steps/loops')
    
    #-----------------------------------------------------------------------
    # LOOP THROUGH THE PERIODS
    #-----------------------------------------------------------------------
    
    for r in range(0,routing_steps):
	     
        print('Rounting Step :',r,' of ', routing_steps)
        
        """"
        --------------------------------------------------------------------
        DEFINE ROUTING STEP INFORMATION 
        -------------------------------------------------------------------
        NOTE 1: days in routing period are variable for qm and monthly
        NOTE 2: for daily routing, EOP date is same as BOP (technically end of day)
        
            ||      DAILY ROUTING       ||       WEEKLY ROUTING
         R  || BOP Date    |  EOP Date  ||   BOP Date   |    EOP Date  
         0  ||  10-1-18       10-1-18   ||    9-28-18         10-4-18
         1  ||  10-2-18       10-2-18   ||    10-5-18        10-11-18
         2  ||  10-3-18       10-3-18   ||   10-12-18        10-18-18
        
            ||      Q-M ROUTING          ||     MONTHLY ROUTING
         R  ||  BOP Date   |  EOP Date   ||  BOP Date   |   EOP Date  
         0  ||   10-1-18      10-8-18    ||   10-1-18       10-31-18
         1  ||   10-9-18      10-15-18   ||   11-1-18       11-30-18
         2  ||  10-16-18      10-23-18   ||   12-1-18       12-31-18       
        
        ITEM OF INTEREST: seconds in routing period (dt)
        """
        if key_intvl=='dy':
            # EOP Date (same day as BOP)
            daysinperiod = 1
            dt = daysinperiod*secinday
            routedateEOP = routedateBOP            

        if key_intvl=='wk':
            # EOP Date (6 days after BOP)
            daysinperiod = 7
            dt = daysinperiod*secinday
            routedateEOP = routedateBOP + db.datetime.timedelta(days=daysinperiod-1)  

        if key_intvl=='qm':
            # EOP Date (6 or 7 days after BOP)
            qtr = db.util.qtr_month_from_date(routedateBOP)
            daysinperiod = db.util.days_in_qtr_mon(routedateBOP.year,routedateBOP.month,qtr)
            dt = daysinperiod*secinday
            routedateEOP = routedateBOP + db.datetime.timedelta(days=daysinperiod-1)

        if key_intvl=='mn':
            # EOP Date (27, 28, 29, 30 days after BOP)
            daysinperiod = db.util.days_in_month(routedateBOP.year,routedateBOP.month)
            dt = daysinperiod*secinday    
            routedateEOP = routedateBOP + db.datetime.timedelta(days=daysinperiod-1)
        
        #--------------------------------------------------------------------
        # ... continue with the calculations below for the routing period ...
        #--------------------------------------------------------------------
        
        info_d = 'BOP: ' + str(routedateBOP) + ' - EOP: '  + str(routedateEOP)
        info_t = ' [dt = ' + str(daysinperiod) + ' d | ' + str(dt) + ' sec]'
        
        # Print out the routing period date information
        print(info_d,info_t)
		 
        # Reset/Initialize the mean level for this period as 1/2 of BOP level.
        spmlvcn = float(spbop * half)
        
        # Reset/Initialize the river flows for this routing step.
        smrrivflow = zero
        
        # Reset/Initialize the 'initial' end of increment level
        speoi = spbop
        
        # Grab the known component data (nts,nbs,ice,div) for the routing step
        spnbs = float(list_nbs_sup[r])
        spdiv = float(list_div_sup[r])
        spice = float(list_icw_sup[r])
        spcsu = float(list_csu_sup[r])
        spgdw = float(list_gdw_sup[r])
    
        # TEMPORARY - TRICK ONLY FOR WHEN FORCED FLOW
        sptemp = float(list_mlv_sup[r])
        smtemp = float(list_flo_smr[r])
        
        # Reset/Initialize the continuity check criteria  
        solveerrorcount = 0
        maxcmserror = solvetestcms
        
        # Calculate number of seconds in an increment for the period
        timesec = dt/spincnumber
        
        # Warn User if potential numeric instability
        if timesec > 222000:
            print('Seconds per Increment is Larger Than Recommended for Period', timesec )
            print('To Resolve: Increase Number of Increments')    
    
        #-----------------------------------------------------------------------
        # LOOP THROUGH THE INCREMENTS IN EACH PERIOD
        #-----------------------------------------------------------------------
 
        for i in range(1,spincnumber+1):
            
            # print('Increment ',i,'out of ',spincnumber)
            
            spboi = speoi
            
            if spsolutionid==1:
                #--------------------------------------------------------------
                #              Solution Method: FORCED (1)
                #--------------------------------------------------------------  
                print('solve forced')
            if spsolutionid==2:
                #--------------------------------------------------------------
                #              Solution Method: PREPROJECT (2)
                #-------------------------------------------------------------- 
                print('solve preproject')
            if spsolutionid==3:
                #--------------------------------------------------------------
                #              Solution Method: PLAN77a (3)
                #--------------------------------------------------------------                 
                print('solve plan77a')
            if spsolutionid==4:
                #--------------------------------------------------------------
                #              Solution Method: PLAN2012 (4)
                #-------------------------------------------------------------- 
                print('solve plan2012')
          
        #----------------------------------------------------------------------
        # OVERVIEW OF VARIABBLES AT END OF LOOP
        #----------------------------------------------------------------------
        
        # At this stage of the routesuperior program...
        # - boi levels are evaluated for final increment in the period
        # - eoi levels are evaluated for final increment in the period
        # - outi flows are mean flows for the final increment in the period
        # - rivflow flows are sum of outi flows over the period
        # - mlv levels are sum of half bop plus each eoi levels
        # - SolveErrors at this stage is # of increments in period when check failed
        
        #-------------------------------------------------------------------------
        # FINAL PERIOD CALCULATIONS
        #-------------------------------------------------------------------------
        # NOT TRUE CALCULATIONS
        # Set EOP values to final EOI level calculated
        speop = 2*sptemp - spbop
       
        # Calculate the final mean lake levels of the period 
        spmlv = sptemp
        
        # Calculate the final mean river flows for the period
        smrrivflow = smtemp
        
        # Calculate the change in levels from bop to eop
        dzsp = speop - spbop
        
        #-------------------------------------------------------------------------
        # CALCULATE THE STORAGE
        #-------------------------------------------------------------------------

        dssp = dzsp*Asp/dt
      
        #-------------------------------------------------------------------------
        # POPULATE LISTS WITH THE RESULTS
        #-------------------------------------------------------------------------
        
        list_mlv_sup[r] = round(spmlv,3)
        list_flo_smr[r] = round(smrrivflow,-2)             
        list_str_sup[r] = round(dssp,-2)
        
        list_dt_bop[r] = routedateBOP
        list_dt_eop[r] = routedateEOP
        list_dt_day[r] = daysinperiod
        list_dt_sec[r] = dt
        #-----------------------------------------------------------------------
        # PREP FOR THE NEXT LOOP/ROUTING PERIOD: RESET DATE,BOP
        #-----------------------------------------------------------------------

        # Reset the new BOP date
        routedateBOP = routedateEOP + oneday
        
        # Reset the BOP levels for the next period as the EOP level from this period
        spbop = speop

    #-----------------------------------------------------------------------
    # POST-CALCULATIONS: SAVE AND DEPOSIT RESULTS INTO VAULT
    #-----------------------------------------------------------------------
    # Deposit ts of mean lake levels (t res is routing period)
    VaultData.deposit_data(kind='mlv', units='m', intvl=key_intvl, loc='sup',
                          first=startdate, last=finaldate, values=list_mlv_sup)
    
    # Deposit ts of mean flows for connecting channels (t res is routing period)
    VaultData.deposit_data(kind='flow', units='cms', intvl=key_intvl, loc='smr', 
                           first=startdate, last=finaldate, values=list_flo_smr)
    #-----------------------------------------------------------------------
    # OUTPUT AND LOG THE CALCULATIONS
    #-----------------------------------------------------------------------
    sp_stor = [list_str_sup]
    sp_file_csv = 'routesuperior_validate.csv'
    
    if spverbosity > 2:
        print('...return 1 storage value, bop, eop, dt_day, dt_sec...')
        sp_log_data = [list_dt_bop,
                       list_dt_eop, 
                       list_dt_day,
                       list_dt_sec,
                       list_str_sup]
        
        sp_name_data = ['bop[date]',
                        'eop[date]',
                        'dt[day]',
                        'dt[sec]',
                        'sup-dst[cms]']    
    if spverbosity > 5:
        print('...return 1 storage value, bop, eop, dt_day, dt_sec...')
        sp_log_data = [list_dt_bop,
                       list_dt_eop, 
                       list_dt_day,
                       list_dt_sec,
                       list_str_sup]
        
        sp_name_data = ['bop[date]',
                        'eop[date]',
                        'dt[day]',
                        'dt[sec]',
                        'sup-mlv[m]',
                        'smr-flw[cms]',
                        'sup-dst[cms]'] 
    
    # TODO: decide what the returned 'value' should be....
    # NOTE - atm this is redundant - since the results are stored in the vault...
    return sp_log_data, sp_name_data, sp_file_csv, sp_stor
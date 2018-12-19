# IMPORTANT: THIS CODE IS NOT COMPLETE. CURRENTLY UNOFFICIAL AND WORK IN PROGRESS!!!!!
#  - sharing with group just so everyone can see how things are evolving
#------------------------------------------------------------------------------
#                             Middle Lake Program
#------------------------------------------------------------------------------
# Coded in python by Zoe Miller zoe.a.miller@usace.army.mil
# Last Edit: 12/19/2018

#------------------------------------------------------------------------------
#                      IMPORT REQUIRED FUNCTIONS
#------------------------------------------------------------------------------

from middle.routemiddle_util import Outflow, GetArea
from middle.solution_iterative import SolveIterative
from middle.solution_matrix import SolveMatrix

import handler.databank as db
   
#----------------------------------------------------------------------------
#                         MAIN PROGRAM: middlelakes
#----------------------------------------------------------------------------

def routemiddle(VaultData, SimulationStartDate, SimulationEndDate, MiddleLakesInfo):
    """
    ----------------------------------------------------------------------
    -                         routemiddle() INPUTS
    ----------------------------------------------------------------------
    Middle Lakes Routing: Solve for period lake levels and flows from the start
    to end date
    
    ABOUT: 
    Function solves the "hydrologic equations" for all periods (routing steps)
    by breaking the problem into incremental steps. 
    
    INPUTS: 4 Total (all required)
    1 object of DataVault type (see handler/databank for more detail)
            - contains time series to be withdrawn
    2 variables of date type:  SimulationStartDate, SimulationEndDate
    1 object of ml_config_info: (see middle_util for more detail)
            - contains crtieria for middlelakes routing and settings
    
    OUTPUT: 4 Lists (at moment - mostly temporary here...to be edited as
    glrrm whole framework is developed)
    
     - ml_log_data : N lists full of routing period 
     - ml_name_data: headers of the ml_log_data
     - ml_file_csv
     - ml_stor

    ----------------------------------------------------------------------
    -                          Program Notes
    ----------------------------------------------------------------------
    PROGRAM HISTORY: 
    
    TRANSLATION: Zoe Miller (USACE Detroit) translated the routemiddle subroutine
    from the CGLRRM into python.  No changes to functionality of program, but
    some additional water balance terms are added (ie consumptive use, 
    groundwater, etc) and the program is written to use the Data Handler.
    - January 2019 - version 1
    
    ORIGINAL SOURCE: Code is translated and adapted from FORTRAN code in the
    CGLRRM. Matt McPherson (formerly USACE Detroit, now HEC) wrote and
    commented the code.  Zoe Miller (USACE Detroit), translated this code to 
    python Original file: midlakes.for and routine name: RouteMiddleLakes. His 
    code was adapted, somewhat, from NOAA GLERL code: 
        https://www.glerl.noaa.gov/pubs/tech_reports/glerl-109/
               
    ----------------------------------------------------------------------
    -                 Overview of Middle Lake Routing Problem
    ----------------------------------------------------------------------    
    WHAT'S THE PROBLEM: "Hyrologic Equations"
    
    How can we solve for future water levels, if the following rates are known 
    variables for each routing time step (period):
        
    - St Marys Flow: QSM (qsm)
    - Net Basin Supplies: NBS (nbs)
    - Ice and Weed Retardation: R (icw) 
    - Diversions: DIV (div)
    - Consumptive Use : CSU (csu)- default is zero
    - Groundwater: GDW (gdw) - default is zero

    ANSWER: Continuity Equations for each Lake
    
    WHAT'S THE BASIC PRINCIPLE: After some period of time dt, The change in
    storage is equal to what enters the lake minus what leaves the lake. 
    Beginning of period is Z0 and end of period is Zt.
 
            dS/dt = Water In - Water Out
        A*(dZ/dt) = Water In - Water Out
   A*(Zt - Z0)/dt = Water In - Water Out 
                0 = Water In - Water Out - A*(Zt - Z0)/dt

    Note that Q = O - R. Connecting Channels are stage-fall minus ice and weed     
    
    mhu:  NBS1 + GDW + QSM - CSU - D1  - Q1 - A1*(dZ1/dt) = 0
    stc:  NBS2 + GDW + Q1  - CSU       - Q2 - A2*(dZ2/dt) = 0
    eri:  NBS3 + GDW + Q2  - CSU - D3  - Q3 - A3*(dZ3/dt) = 0
    
    HOW'S THE PROBLEM SOLVED: Numerical Methods
    
    ----------------------------------------------------------------------
    -                       Middle Lakes Routing Steps
    ----------------------------------------------------------------------
    LINGO:  BOP = Beginning of Period,
            EOP = End of Period,
            BOI = Beginning of Increment,
            EOI = End of of Increment, 
            
    Depending on how the equations are rearranged, they are of a format, where
    either dZ's or Zt's can be solved in at least two different ways for some 
    interval.  In this problem, break up each period into smaller increments.
    For each increment, solve hydrologic problem, to get EOI levels and flows.

    CURRENTLY TWO SOLUTION METHODS: Iterative and Matrix
 
    STEP 1: Calculations for one period at a time (loop through routing periods)
         - Knowns from start or prior period: BOP, solvererrorcount,
    STEP 2: Calculations for one increment at a time (loop through increments)
         - Knowns: BOI, MLVCN
         - Solve for EOI levels using a solution method
         - Caclulate EOI outflows
         - Calculate cumulative sum of EOI levels and connecting channel flows
             purely a coding tactic for averaging later for routing period mlv
         - Evaluate performance of solution with conservation of mass check
    STEP 3: Final Calculations for a period, after incremental calculations
         - Set EOP equal to final EOI
         - Calculate mean levels for routing period
         - Reset any variables for next routing period
    
    TIME GRID EXAMPLE:
 
    | inc1   inc2   inc3   inc4   inc5 | inc1  inc2   inc3   inc4   inc5 | ...
    *------*------*------*------*------*-----*------*------*------*------* ...
    |                                  |                                 | ...
    *----------------------------------*---------------------------------*
                 period1                            period2
     
    Variables of interest are the routing period mean levels and flows
 
    Key Averaging Trick: Averaging the mean levels of the increments, over the 
    period can be rearranged to the following:
        
             MLV = [0.5*BOP + sumEOI + 0.5*EOP]/n
             
    NOTE:  sumEOI is only for 1 through n-1 whre n is number if increments !
    
    ----------------------------------------------------------------------
    -                     Extra Sources & Reading
    ----------------------------------------------------------------------
    - CGLRRM_Update.docx (Zoe ...working on documentation...)
    - Quinn, F.H. Hydrologic response model of the North American Great Lakes. 
          J. of Hydrology 37:295-307 (1978).
    -  CGLRRM draft manual (February 2001)   
    """
    #-----------------------------------------------------------------------
    # CHECK INPUTS ENTERED INTO THE PROGRAM
    #-----------------------------------------------------------------------
    # TODO: add a check outside of midlakes for data existing for sim length
    if not VaultData:
        raise Exception('Middle Lakes Program missing VaultData')
    if not SimulationStartDate:
        raise Exception('Middle Lakes Program missing SimulationStartDate')
    if not SimulationEndDate:
        raise Exception('Middle Lakes Program missing SimulationEndDate')
    if not MiddleLakesInfo:
        raise Exception('Middle Lakes Program missing MiddleLakesInfo')
    
    #-----------------------------------------------------------------------
    # UNPACK VARIABLES NEEDED FOR SOLUTION CALCULATIONS
    #-----------------------------------------------------------------------
    
    # "Unpack" the 2 dates of formate (date time)
    startdate = SimulationStartDate
    finaldate = SimulationEndDate

    # define the constants (float)
    zero = 0.0
    half = 0.5
    secinday = 86400
  
    # Unpack from the ml_config_info
    
    # period info
    key_intvl = MiddleLakesInfo.keyintvl         # str for DataVault retrieval
    oneday = MiddleLakesInfo.timedelta_one_day                # type timedelta
    
    # increment info
    mlincnumber = MiddleLakesInfo.mlincrements                      # integer 
    mlsolutionid = MiddleLakesInfo.solutionindicator                # integer
    
    # for logging (integer -- logging and handling to be added)
    mlverbosity = MiddleLakesInfo.mlverbosity                        

    # for continuity check (todo: what value should this really be?)
    solvetestcms = 75	
    
    # rounging values (integer, potentially negative)
    roundflo = MiddleLakesInfo.roundflo
    roundlev = MiddleLakesInfo.roundlev

    # stage-fall discharge coefficients (floats)
    mhk  = MiddleLakesInfo.stagefall_mhk
    mhym = MiddleLakesInfo.stagefall_mhym
    mha  = MiddleLakesInfo.stagefall_mha
    mhb  = MiddleLakesInfo.stagefall_mhb
    mhwt = MiddleLakesInfo.stagefall_mhwt
    mhc  = MiddleLakesInfo.stagefall_mhc
	
    sck  = MiddleLakesInfo.stagefall_sck
    scym = MiddleLakesInfo.stagefall_scym
    sca  = MiddleLakesInfo.stagefall_sca
    scb  = MiddleLakesInfo.stagefall_scb
    scwt = MiddleLakesInfo.stagefall_scwt
    scc  = MiddleLakesInfo.stagefall_scc
	
    erk  = MiddleLakesInfo.stagefall_erk
    erym = MiddleLakesInfo.stagefall_erym
    era  = MiddleLakesInfo.stagefall_era
    erb  = MiddleLakesInfo.stagefall_erb
    erwt = MiddleLakesInfo.stagefall_erwt
    erc  = MiddleLakesInfo.stagefall_erc
     
    # starting levels (floats)
    start_lev_mhu = MiddleLakesInfo.startmhu
    start_lev_stc = MiddleLakesInfo.startstc
    start_lev_eri = MiddleLakesInfo.starteri
 
    # Unpack from the DataVault [requires: StartDate,EndDate & key_intvl]
    
    try:
        # Withdraw St. Mary's River [list]
        list_flo_smr = VaultData.withdraw(kind='flow', units='cms',
            intvl=key_intvl, loc='smr', first=startdate, last=finaldate).dataVals
         
        # define the number of routing period based on St. Mary's list (integer)
        list_len = len(list_flo_smr)
    except:
        raise Exception('Middle Lakes: Missing St Marys Flow at Interval =',key_intvl)
    
    try:
        # Withdraw Net Basin Supplies  [list] 
        list_nbs_mhu = VaultData.withdraw(kind='nbs', units='cms', 
            intvl=key_intvl, loc='mhu', first=startdate, last=finaldate).dataVals
        list_nbs_stc = VaultData.withdraw(kind='nbs', units='cms',
            intvl=key_intvl, loc='stc', first=startdate, last=finaldate).dataVals
        list_nbs_eri = VaultData.withdraw(kind='nbs', units='cms',
            intvl=key_intvl, loc='eri', first=startdate, last=finaldate).dataVals
    except:
        raise Exception('Middle Lakes: Missing Net Basin Supplies')
    
    try:
        # Withdraw Ice and Weed Retardation [list] 
        list_icw_mhu = VaultData.withdraw(kind='icw', units='cms',
            intvl=key_intvl, loc='scr', first=startdate, last=finaldate).dataVals
        list_icw_stc = VaultData.withdraw(kind='icw', units='cms',
            intvl=key_intvl, loc='det', first=startdate, last=finaldate).dataVals
        list_icw_eri = VaultData.withdraw(kind='icw', units='cms',
            intvl=key_intvl, loc='nia', first=startdate, last=finaldate).dataVals
    except:
        raise Exception('Middle Lakes: Missing Ice and Weed Retardation Values')
   
    try:
        # Withdraw Diversions [list] 
        list_div_mhu = VaultData.withdraw(kind='flow', units='cms', 
            intvl=key_intvl, loc='chi', first=startdate, last=finaldate).dataVals
        list_div_eri = VaultData.withdraw(kind='flow', units='cms',
            intvl=key_intvl, loc='wel', first=startdate, last=finaldate).dataVals
    except:
        # Create Diverion lists of constants if none in Vault [list]
        # TODO: Add to config file
        list_div_mhu = [91.0]*list_len
        list_div_eri = [190.0]*list_len
        print('Middle Lakes: Diversions Missing from Vault, set to 91.0 & 190.0')
  
    try:
        # Withdraw Consumptive Use [list]
        list_csu_mhu = VaultData.withdraw(kind='con', units='cms',
            intvl=key_intvl, loc='mhu', first=startdate, last=finaldate).dataVals
        list_csu_stc = VaultData.withdraw(kind='con', units='cms',
            intvl=key_intvl, loc='stc', first=startdate, last=finaldate).dataVals
        list_csu_eri = VaultData.withdraw(kind='con', units='cms', 
            intvl=key_intvl, loc='eri',  first=startdate, last=finaldate).dataVals        
    
    except:
        # Create Consumptive Use lists of constants if none in Vault [list]
        list_csu_mhu = [0]*list_len
        list_csu_stc = [0]*list_len
        list_csu_eri = [0]*list_len
        print('Middle Lakes: Consumptive Use Missing from Vault, set to zeros')
    
      
    try:
        # Withdraw Groundwater   
        list_gdw_mhu = VaultData.withdraw(kind='flw', units='cms',
            intvl=key_intvl, loc='mhu', first=startdate, last=finaldate).dataVals
        list_gdw_stc = VaultData.withdraw(kind='flw', units='cms',
            intvl=key_intvl, loc='stc', first=startdate, last=finaldate).dataVals
        list_gdw_eri = VaultData.withdraw(kind='flw', units='cms', 
            intvl=key_intvl, loc='eri', first=startdate, last=finaldate).dataVals        
    
    except:
        # Create Groundwater lists of constants if none in Vault [list]
        list_gdw_mhu = [0]*list_len
        list_gdw_stc = [0]*list_len
        list_gdw_eri = [0]*list_len
        print('Middle Lakes: Groundwater Missing from Vault, set to zeros')   
   
    try:
        # Withdraw CGIP levels if provided (else devault value: 171.16 )
        list_mlv_gip = VaultData.withdraw(kind='mlv', units='m',
            intvl=key_intvl, loc='nia', first=startdate, last=finaldate).dataVals
    except:
        # Create CGIP levels of constants if none in Vault [list]
        list_mlv_gip = [171.16]*list_len  
        print('Middle Lakes: CGIP Missing from Vault, set to 171.16 constant') 

    # COMPLETED: UNPACKING VARIABLES NEEDED FOR SOLUTION CALCULATIONS
    
    #-----------------------------------------------------------------------
    # BEGIN MIDDLE LAKES ROUTING PROGRAM
    #-----------------------------------------------------------------------

    #-----------------------------------------------------------------------
    # INITIALIZATION
    #-----------------------------------------------------------------------
    print('The Start Date: ', startdate)
    print('The End Date: ', finaldate)

    # Calculate the net total supply for michigan-huron & create empty lists
    
    list_nts_mhu = [None]*list_len
    
    for s in range(0,list_len):
        list_nts_mhu[s] = list_flo_smr[s] + list_nbs_mhu[s]
        
    list_mlv_mhu = [None]*list_len
    list_mlv_stc = [None]*list_len
    list_mlv_eri = [None]*list_len
	
    list_flo_src = [None]*list_len
    list_flo_det = [None]*list_len
    list_flo_nia = [None]*list_len
	 
    list_str_mhu = [None]*list_len
    list_str_stc = [None]*list_len
    list_str_eri = [None]*list_len 
    
    list_dt_bop = [None]*list_len
    list_dt_eop = [None]*list_len
    list_dt_day = [None]*list_len
    list_dt_sec = [None]*list_len
    
    # Initialize BOP levels as starting levels for first period
    mhbop = start_lev_mhu
    scbop = start_lev_stc
    erbop = start_lev_eri
    cgbop = 171.16 # pg 24 of cglrrm draft manual ... or  list_mlv_gip[0] ?
	
    # Initiate the routing date as the start date for the first period
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
        mhmlvcn = float(mhbop * half)
        scmlvcn = float(scbop * half)
        ermlvcn = float(erbop * half)
        
        # Reset/Initialize the river flows for this routing step.
        stcrivflow = zero
        detrivflow = zero
        niarivflow = zero
        
        # Reset/Initialize the 'initial' end of increment level
        mheoi = mhbop
        sceoi = scbop
        ereoi = erbop
        cgeoi = cgbop                    # just needed for nia equations
        
        # Grab the known component data (nts,nbs,ice,div) for the routing step
        mhnts = float(list_nts_mhu[r] - list_div_mhu[r])
        scnbs = float(list_nbs_stc[r])
        ernbs = float(list_nbs_eri[r]- list_div_eri[r])
        
        mhice = float(list_icw_mhu[r])
        scice = float(list_icw_stc[r])
        erice = float(list_icw_eri[r])
        
        mhdiv = float(list_div_mhu[r])
        erdiv = float(list_div_eri[r])
        
        mhcsu = float(list_csu_mhu[r])
        sccsu = float(list_csu_stc[r])
        ercsu = float(list_csu_eri[r])
     
        mhgdw = float(list_gdw_mhu[r])
        scgdw = float(list_gdw_stc[r])
        ergdw = float(list_gdw_eri[r])
        
        # Reset/Initialize the continuity check criteria  
        solveerrorcount = 0
        maxcmserror = solvetestcms
        
        # Calculate number of seconds in an increment for the period
        timesec = dt/mlincnumber
        
        # Warn User if potential numeric instability
        if timesec > 222000:
            print('Seconds per Increment is Larger Than Recommended for Period', timesec)
            print('To Fix: Increase Number of Increments')
        
        #-----------------------------------------------------------------------
        # LOOP THROUGH THE INCREMENTS IN EACH PERIOD
        #-----------------------------------------------------------------------
 
        for i in range(1,mlincnumber+1):
            
            print('Increment ',i,'out of ',mlincnumber)
            
            mhboi = mheoi
            scboi = sceoi
            erboi = ereoi
            cgboi = cgeoi
            
            if mlsolutionid==1:
                #--------------------------------------------------------------
                #              Solution Method: ITERATIVE (1)
                #--------------------------------------------------------------            
                
                [mheoi,sceoi,ereoi] = SolveIterative(BOI_MHU = mhboi,
                                                     BOI_STC = scboi,
                                                     BOI_ERI = erboi,
                                                     BOI_CGI = cgboi,
                                                     NTS_MHU = mhnts,
                                                     NBS_STC = scnbs,
                                                     NBS_ERI = ernbs,
                                                     ICW_MHU = mhice,
                                                     ICW_STC = scice,
                                                     ICW_ERI = erice,
                                                     DIV_MHU = mhdiv,
                                                     DIV_ERI = erdiv,
                                                     deltat = timesec,
                                                     MiddleLakesInfo = MiddleLakesInfo)
                cgeoi = cgboi
            
            elif mlsolutionid==2:
                #--------------------------------------------------------------
                #               Solution Method: MATRIX (2)
                #--------------------------------------------------------------
            
                [mheoi,sceoi,ereoi] = SolveMatrix(BOI_MHU = mhboi,
                                                     BOI_STC = scboi,
                                                     BOI_ERI = erboi,
                                                     BOI_CGI = cgboi,
                                                     NTS_MHU = mhnts,
                                                     NBS_STC = scnbs,
                                                     NBS_ERI = ernbs,
                                                     ICW_MHU = mhice,
                                                     ICW_STC = scice,
                                                     ICW_ERI = erice,
                                                     DIV_MHU = mhdiv,
                                                     DIV_ERI = erdiv,
                                                     deltat = timesec,
                                                     MiddleLakesInfo = MiddleLakesInfo)
                cgeoi = cgboi
            
            #----------------------------------------------------------------------
            # CALCULATE CONNECTING CHANNEL FLOW FOR THE INCREMENT
            #----------------------------------------------------------------------
            """
            BOP                       MLV                          EOP           
             |                       period1                        |                                                        
             *------------------------------------------------------*-- ...            
             |       inc1             inc2               inc3       |    
             *-----------------*-----------------*------------------*-- ....
             BOI1             EOI1
                              BOI2              EOI2
                                                BOI3              EOI3
            
            How is the connecting channel flow calculated for the increment?
                - stage-fall discharge equations with ice and weed retardation
                - evaluate at the mean increment levels: Zbar = 0.5(boi + eoi)
            	
            NOTE: Niagara river flow (erouti) has cgmlv as ZD, which would be the mean
            of the downstream body of water, in this case the chippewa-grass island
            pool (cgmlv).  A value is required for ZD to solve for a two gage equation,
            however while a single gage equation is used, an input is still required, 
            even though any value can technically be entered without influencing
            the result, because of the fact that ZD^0 = 1..
            """
            
            # Calculate mean level for the increment
            mhmlvi = half*(mhboi + mheoi)
            scmlvi = half*(scboi + sceoi)
            ermlvi = half*(erboi + ereoi)
            cgmlvi = half*(cgboi + cgeoi)
            
            # Mean Connecting Channel Flow for the Increment
            mhouti = Outflow(mhmlvi,scmlvi,mhk,mhym,mha,mhb,mhwt,mhc) - mhice
            scouti = Outflow(scmlvi,ermlvi,sck,scym,sca,scb,scwt,scc) - scice
            erouti = Outflow(ermlvi,cgmlvi,erk,erym,era,erb,erwt,erc) - erice

            # Accumulate flows over all increments in a period (for later averaging)  
            stcrivflow = stcrivflow + mhouti
            detrivflow = detrivflow + scouti
            niarivflow = niarivflow + erouti
            
            #----------------------------------------------------------------------
            # HOW WELL DID SOLUTION SOLVE THE PROBLEM: conservation of mass
            #----------------------------------------------------------------------
            """
            Recall that
                        Water In =  Water Out + A*(dz/dt)
                               0 =  Water Out - Water In + A*(dZ/dt)
                               0 =  Water Out - Water In + A*((EOI - BOI)/dt)
            
            We used a numerical method to solve the problem (for EOI and
            therefore dZ), so it probably won't be an exact solution. 
            
            How well did we do? [in cms]
            To check, substitue in the values calcualted (EOI), and see if the
            RHS of the equation is ~ 0. In a perfect world, this would be = 0.
            """
            # calculate the difference (delta Z) in meters from eoi and boi
            dzmhi = mheoi - mhboi
            dzsci = sceoi - scboi
            dzeri = ereoi - erboi
            
            # calculate the surface areas for each lake (generally constant)
            Amh = GetArea('mh',mheoi)
            Asc = GetArea('sc',sceoi)
            Aer = GetArea('er',ereoi)
            
            # calculate the change in storage dS/dt = A*dz/dt (for increment)
            dsmhi = dzmhi*Amh/timesec
            dssci = dzsci*Asc/timesec
            dseri = dzeri*Aer/timesec
            
            # need to add in the consumptive use and groundwater terms
            
            # calculate the mass balance (cms): O - I + ds/dt
            mhsolutioncheck = (mhouti + mhdiv + mhcsu) - (mhnts          + mhgdw) + dsmhi
            scsolutioncheck = (scouti         + sccsu) - (scnbs + mhouti + scgdw) + dssci
            ersolutioncheck = (erouti + erdiv + ercsu) - (ernbs + scouti + ergdw) + dseri
            
            # mhsolutioncheck = mhouti          - mhnts + mhdiv + dsmhi
            # scsolutioncheck = scouti - mhouti - scnbs         + dssci
            # ersolutioncheck = erouti - scouti - ernbs + erdiv + dseri
            
            mhsolutioncheck = round(abs(mhsolutioncheck),roundflo)
            scsolutioncheck = round(abs(scsolutioncheck),roundflo)
            ersolutioncheck = round(abs(ersolutioncheck),roundflo)
            
            #------------------------------------------------------------------
            # DOES THE SOLUTION MEET CONSERVATION OF MASS STANDARDS?
            #------------------------------------------------------------------
            """ 
            User/modeler defined what threshold is appropirate for meeting 
            conservation of mass [cms].  What values (100 cms, 10 cms) will
            be close enough to 0?
            
            TODO: Zoe needs to better understand how these thresholds are 
            defined -- and make sure she is using an appropirate threshold
            for the increment (here) and ultimately routing period (below)!
            """
            
            # Recall the user defined solvetestcms (75). 
            solvetestcms

            if ersolutioncheck >= solvetestcms or scsolutioncheck >= solvetestcms or mhsolutioncheck >= solvetestcms:
                # Count the number of increments in the period where conservation
                # of mass check fails for at least one connecting channel
                solveerrorcount = solveerrorcount + 1
               
                # Track largest magnitude (cms) of when cons. of mass fails
                if ersolutioncheck > maxcmserror:
                    # New continuity error (cms) record for period
                    maxcmserror = ersolutioncheck
                    print('RECORD ERROR [PRODUCED BY ER]: ',maxcmserror)
                if scsolutioncheck > maxcmserror:
                    # New continuity error (cms) record for period
                    maxcmserror = scsolutioncheck
                    print('RECORD ERROR [PRODUCED BY SC]: ',maxcmserror)
                if mhsolutioncheck > maxcmserror:
                    # New continuity error (cms) record for period
                    maxcmserror = mhsolutioncheck
                    print('RECORD ERROR [PRODUCED BY MH]: ',maxcmserror)
            else:
                # SCENARIO: conservation of mass check PASSES for an increment
                print('FOR THIS INCREMENT[',i,']IN PERIOD[',r,'] ALL CONS. OF MASS CHECKS PASS')
                # In other words, each lake meets the following criteria
                #    	O - I + aveDELTAZ/dt < SolutionTest
 
            #------------------------------------------------------------------
            # CALCULATE CUMULATIVE MEAN LEVEL AT EACH INCREMENT
            #------------------------------------------------------------------
            """
            What is MLVCN? It is the cumulative sum of the 0.5* BOP plus each EOI
            Why calculate? It will be the numerator for the averaging at EOP
            
            NOTE that prior to the first increment MLV was set to half BOP
            NOTE here it is clear that each EOI is added to the cumulative sum
            
            NOTE that in the final increment, this means the final EOI 
            is added. This means in the final increment we added double the 
            amount supposed to be added (0.5*EOP=0.5*EOI). After the loop
            exists, one half of the EOP will be removed. 
            """
            mhmlvcn = mhmlvcn + mheoi
            scmlvcn = scmlvcn + sceoi
            ermlvcn = ermlvcn + ereoi
                 
        #----------------------------------------------------------------------
        # OVERVIEW OF VARIABLES AT END OF LOOPING THROUGH INCREMENTS
        #----------------------------------------------------------------------

        # - boi levels are evaluated for final increment in the period
        # - eoi levels are evaluated for final increment in the period
        # - outi flows are mean flows for the final increment in the period
        # - rivflow flows are sum of outi flows over the period
        # - mlv levels are sum of half bop plus each eoi levels
        # - SolveErrors at this stage is # of increments in period when check failed
        
        #-------------------------------------------------------------------------
        # FINAL PERIOD CALCULATIONS: LEVEL, FLOWS & STORAGE
        #-------------------------------------------------------------------------
            
        # Set EOP values to final EOI level calculated
        mheop = mheoi
        sceop = sceoi
        ereop = ereoi
        
        """
        BOP                      MLV                          EOP           
        |                       period1                        |                                                        
        *------------------------------------------------------*-- ...            
        |       inc1             inc2               inc3       |    
        *-----------------*-----------------*------------------*-- ....
              MLV_I1             MLV_I2             MLV_I3
          0.5(BOI1+ EOI1)    0.5(BOI2+EOI2)     0.5(BOI2+EOI2)
                                                 
        
        MEAN OF THE PERIOD = AVERAGE OF MEAN INCREMENTAL VALUES
        NUMERATOR = "Sum of Mean Increment Levels"
        NUMERATOR = 0.5*BOP + sumEOI + 0.5*BOP     [EOI only from 1 to n-1]!
        DENOMINATOR = n
        MLV = [0.5BOP + sumEOI - 0.5BOP]/n      [when EOI from 1 through n]
        """
        # Ensure the numerator is the correct value: substract extra half of eop
        mhmlvcn = mhmlvcn - mheop * half
        scmlvcn = scmlvcn - sceop * half
        ermlvcn = ermlvcn - ereop * half
        
        # Calculate the final mean lake levels of the period 
        mhmlv = mhmlvcn/mlincnumber
        scmlv = scmlvcn/mlincnumber
        ermlv = ermlvcn/mlincnumber
        
        # Calculate the final mean river flows for the period
        stcrivflow = stcrivflow/mlincnumber
        detrivflow = detrivflow/mlincnumber
        niarivflow = niarivflow/mlincnumber
        
        # Calculate the change in levels from bop to eop
        dzmh = mheop - mhbop
        dzsc = sceop - scbop
        dzer = ereop - erbop
        
        #-------------------------------------------------------------------------
        # CALCULATE THE STORAGE
        #-------------------------------------------------------------------------

        dsmh = dzmh*Amh/dt
        dssc = dzsc*Asc/dt
        dser = dzer*Aer/dt
      
        #-------------------------------------------------------------------------
        # POPULATE LISTS WITH THE RESULTS
        #-------------------------------------------------------------------------
        
        list_mlv_mhu[r] = round(mhmlv,roundlev)
        list_mlv_stc[r] = round(scmlv,roundlev)
        list_mlv_eri[r] = round(ermlv,roundlev)

        list_flo_src[r] = round(stcrivflow,roundflo)
        list_flo_det[r] = round(detrivflow,roundflo)     
        list_flo_nia[r] = round(niarivflow,roundflo)
        
        list_str_mhu[r] = round(dsmh,roundflo)
        list_str_stc[r] = round(dssc,roundflo)
        list_str_eri[r] = round(dser,roundflo)
        
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
        mhbop = mheop
        scbop = sceop
        erbop = ereop

    #-----------------------------------------------------------------------
    # POST-CALCULATIONS: SAVE AND DEPOSIT RESULTS INTO VAULT
    #-----------------------------------------------------------------------
    
    # Deposit ts of mean lake levels (t res is routing period)
    VaultData.deposit_data(kind='mlv', units='m', intvl=key_intvl, loc='mhu',
                          first=startdate, last=finaldate, values=list_mlv_mhu)
        
    VaultData.deposit_data(kind='mlv', units='m', intvl=key_intvl, loc='stc', 
                           first=startdate, last=finaldate, values=list_mlv_stc)
    
    VaultData.deposit_data(kind='mlv', units='m', intvl=key_intvl, loc='eri', 
                           first=startdate, last=finaldate, values=list_mlv_eri)
	
    # Deposit ts of mean flows for connecting channels (t res is routing period)
    VaultData.deposit_data(kind='flow', units='cms', intvl=key_intvl, loc='scr', 
                           first=startdate, last=finaldate, values=list_flo_src)
    
    VaultData.deposit_data(kind='flow', units='cms', intvl=key_intvl, loc='det', 
                           first=startdate, last=finaldate, values=list_flo_det)
    
    VaultData.deposit_data(kind='flow', units='cms', intvl=key_intvl, loc='nia', 
                           first=startdate, last=finaldate, values=list_flo_nia)

    # Deposit ts of change in storage (t res is routing period)
    #VaultData.deposit_data(kind='flow', units='cms', intvl=key_intvl, loc='mhu', 
    #                       first=startdate, last=finaldate, values=list_str_mhu)
    
    #VaultData.deposit_data(kind='flow', units='cms', intvl=key_intvl, loc='stc', 
    #                       first=startdate, last=finaldate, values=list_str_stc)
    
    #VaultData.deposit_data(kind='flow', units='cms', intvl=key_intvl, loc='eri', 
    #                       first=startdate, last=finaldate, values=list_str_eri)
    
    #-----------------------------------------------------------------------
    # OUTPUT AND LOG THE CALCULATIONS
    #-----------------------------------------------------------------------
    ml_stor = [list_str_mhu, list_str_stc, list_str_eri]
    ml_file_csv = 'routemiddle_validate.csv'
    
    if mlverbosity > 2:
        print('...return 3 storage values, bop, eop, dt_day, dt_sec...')
        ml_log_data = [list_dt_bop,
                       list_dt_eop, 
                       list_dt_day,
                       list_dt_sec,
                       list_str_mhu, list_str_stc,list_str_eri]
        
        ml_name_data = ['bop[date]',
                        'eop[date]',
                        'dt[day]',
                        'dt[sec]',
                        'mhu-dst[cms]','stc-dst[cms]','eri-dst[cms]']
        
    if mlverbosity > 5:
        print('...return 3 levels, 3 flows, 3 storage values, bop, eop, dt_day, dt_sec...')
        ml_log_data = [list_dt_bop,
                       list_dt_eop, 
                       list_dt_day,
                       list_dt_sec,
                       list_mlv_mhu,list_mlv_stc,list_mlv_eri,
                       list_flo_src,list_flo_det,list_flo_nia,
                       list_str_mhu,list_str_stc,list_str_eri]
        
        ml_name_data = ['bop[date]',
                        'eop[date]',
                        'dt[day]',
                        'dt[sec]',
                        'mhu-mlv[m]','stc-mlv[m]','eri-mlv[m]',
                        'scr-flw[cms]','det-flw[cms]','nia-flw[cms]',
                        'mhu-dst[cms]','stc-dst[cms]','eri-dst[cms]']

    return ml_log_data, ml_name_data, ml_file_csv, ml_stor

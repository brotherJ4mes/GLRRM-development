# IMPORTANT: THIS CODE IS INCOMPLETE AND UNOFFICIAL AND WORK IN PROGRESS!!!!!
#  - sharing with group just so everyone can see how things are evolving
#------------------------------------------------------------------------------
#                          Middle Lake Solution: Iterative
#------------------------------------------------------------------------------
# Coded in python by Zoe Miller zoe.a.miller@usace.army.mil
# Last Edit: 12/19/2018
#
#------------------------------------------------------------------------------
#                  What's the Problem: "Hydrologic Equations"
#------------------------------------------------------------------------------
# Note that Q = O - R. Connecting Channels are stage-fall minus ice and weed 
#
# MASTER SET OF EQUATIONS: HYDROLOGIC MODEL
# mhu:  NBS1 + QSM - D1  - Q1 - A1*(dZ1/dt) = 0
# stc:  NBS2 + Q1        - Q2 - A2*(dZ2/dt) = 0
# eri:  NBS3 + Q2  - D3  - Q3 - A3*(dZ3/dt) = 0
# 
# PROBLEM: NEED TO SOLVE FOR UNKNOWNS Z(t=1):
#                      Zmhu,Zstc,Zeri at t=1 and/or 
#                      dZmhu,dZstc,dZeri for dt
#
# WHEN: THE INITIAL VALUES OF Z(t=0) ARE KNOWN:
#                      Zmhu,Zstc,Zeri are known at t=0. [Initial Values]
#
# ---------------------------------------------------------------------------
#            How to Solve the Problem: Iterative Solution
# ---------------------------------------------------------------------------
# ABOUT ITERATIVE SOLUTION:  

# 1) Recall that NBS, D, R & QO's are rates over the time step (dt)
# 2) Observe that when the true solution (eoi values) are known for the system 
#        of equations then the equations above holds true (equal zero). 
#        Conversevely, when a set of values are pluged into those equations 
#        that are not the true solution, the equations will not equal zero.
#        
#        Background: think of the simple problem of solving for the square
#        root of 3.  What value z solves for f(z) = 3^(0.5) - z = 0 ? 
#        The true solution is 1.732... , therefore 3^(0.5) - 1.732 = 0 (close enough)
#        A wrong solution is 2 , therefore 3^(0.5) - 2 != 0 (it equals -.2679)
# 3) Note that the system of equations above hold true when the LHS is evaluated
#        over the period. NBS,R,D are assumed constant over a period. The best
#        approximation of Q over the period can be calcuated using the mean
#        level's over the period: Q(meanZup,meanZdown)
# 4) About Newton-Raphson or Newton's method: This is a guessing game
# 5) Select a set of values to make a first guess (first iteration) for eoi
#           - calculate what mlv would be for such an eoi
#           - calculate what the flows would be such an eoi
#           - check how well the guess was ... how do we do that?
# 6) Check how well the guess was:
#        An important rule in Newton's method is that when the true solution
#        is guessed, then the Q's can be evaluated with the mean Z's and then
#        all the terms on the LHS above are known (KV). Therefore
#              KV = A*(dZ/dt) = A *(zeoi - zboi)/dt
#              KV*(dt/A) = zeoi - zboi 
#              zeoi = KV*(dt/A) + zboi
#         
#        Zmhu,eoi = [NTSmhu - Dchi - Qscr       ]*(dt/Amhu) + Zmhu,boi
#        Zstc,eoi = [NBSstc + Qscr - Qdet       ]*(dt/Astc) + Zstc,boi
#        Zeri,eoi = [NBSeri - Dwel + Qdet - Qnia]*(dt/Aeri) + Zeri,boi 
# 
# 7) How did we check how well the guess was? Well, technically we calculated
#        a new eoi value. But if we had guessed the correct one, then we 
#        just would have returned the value we guessed. Make sense?
# 8) Now let's quantify how much error there is between guesses (iterations).
#        Recall, when the true solution is guessed, this in theory that the error
#        from iteration to iteration would be zero.
#        ERROR = |GUESS - LAST| = |TRUE - LAST| = |TRUE - TRUE| = 0 
# 9) OK, well, technically we are solving for the eoi levels (z). However,
#        the flows are dependent on the mean level for the increment. Due
#        to some reasoning (ZM is hypothesizing here...validate later), likely
#        that the flow values are much more variable than the levels and of 
#        a greater magnitude, we will set our error tolerance for the flows
#        and not the actual levels we solved for.
# 10) So, that means recalculating the flows, at the mean level for the increment
#        and comparing that flow to the flow from the previous iteration.
# 11) Once, the errors are within a certain threshold for each outflow, then 
#        it is safe to assume that we have converged upon the solution.
#             when        
#        ERROR = |GUESS - LAST|  <  TOLERANCE
#             accept GUESS as TRUE SOLUTION
# 12) Finally, one more step is added here, to address when the solution is 
#        oscilating around the true solution, but not meeting the flow criteria
# 13) Once the solution's error ('diff') is met for convergence criteria, then
#        there is no need to keep iterating - that would be a waste of time
#        so stop iterating and accept the latest eoi levels as the true solution  



# TODO ZM needs to learn about Matt's criteria for these thresholds.
# Seems like there should be guidance on what is appropriate for the period
# and ultimately increment. Recall user defined thresholds from config file.
# TODO ZM doesn't get how -999.99 are appropriate first guesses... ZM not getting 
# this to work always - probably becuase of influenc eof outflow equations and
# the checks that exist in there...
# -----------------------------------------------------------------------------
#                 Import python modules/packages or functions
#------------------------------------------------------------------------------

from middle.routemiddle_util import Outflow, GetArea 

#----------------------------------------------------------------------------
#                         Function: SolveIterative
#----------------------------------------------------------------------------
        
def SolveIterative(BOI_MHU=None,BOI_STC=None,BOI_ERI=None,BOI_CGI=None,
                NTS_MHU=None,NBS_STC=None,NBS_ERI=None,
                ICW_MHU=None,ICW_STC=None,ICW_ERI=None,
                DIV_MHU=None             ,DIV_ERI=None,
                CSU_MHU=None,CSU_STC=None,CSU_ERI=None,
                GDW_MHU=None,GDW_STC=None,GDW_ERI=None,
                deltat = None,
                MiddleLakesInfo=None):
    """
    ----------------------------------------------------------------------
    -                     SolveIterative() INPUTS
    ----------------------------------------------------------------------

    Middle Lakes Solution: Iterative Solution
    
    ABOUT: 
    Function solves the "hydrologic equations" for one increment interval using
    the iterative solution.
    
    INPUTS:
    12 float variables: 
    boi: mhu, stc, eri, cgi
    nts: mhu
    nbs: stc, eri
    icw: mhu, stc, eri
    div: mhu, eri
    & 1 variable object of class/type "ml_config_inf". see middle_util.py for more details.
    MiddleLakesInfo
    
    ----------------------------------------------------------------------
    -                          Program Notes
    ----------------------------------------------------------------------
    PROGRAM HISTORY: 
    
    TRANSLATION: Zoe Miller (USACE Detroit) translated the IterativeRoute subroutine
    from the CGLRRM into python.  No changes to functionality of program, but
    some additional water balance terms are added (ie consumptive use, 
    groundwater, etc).
    - January 2019 - version 1
    
    ORIGINAL SOURCE: Code is translated and adapted from FORTRAN code in the
    CGLRRM. Matt McPherson (formerly USACE Detroit, now HEC) wrote and
    commented the code.  Zoe Miller (USACE Detroit), translated this code to 
    python.
               
    ----------------------------------------------------------------------
    -                Overview of Middle Lake Routing Problem
    ----------------------------------------------------------------------    
    WHAT'S THE PROBLEM: "Hyrologic Equations"
    
    How can we solve for future water levels, if the following rates are known 
    variables for each routing time step (increment):
        
    - St Marys Flow: QSM (qsm)
    - Net Basin Supplies: NBS (nbs)
    - Ice and Weed Retardation: R (icw) 
    - Diversions: DIV (div)
    - Consumptive Use : CSU (csu)- default is zero
    - Groundwater: GDW (gdw) - default is zero
    
    - mhu BOI level 
    - stc BOI level
    - eri BOI level

    KEY ASSUMPTION MADE TO EMPLOY ITERATIVE SOLUTION:
        flow rate for interval assumed to be equal to flow evaluated at mean
        of boi and eoi levels
    
    
        
    
    ----------------------------------------------------------------------
    -                       Iterative Solution Steps
    ----------------------------------------------------------------------    
    


    ----------------------------------------------------------------------
    -                     Extra Sources & Reading
    ----------------------------------------------------------------------
    - CGLRRM_Update.docx (Zoe ...working on documentation...)
    - Quinn, F.H. Hydrologic response model of the North American Great Lakes. 
          J. of Hydrology 37:295-307 (1978).
    -  CGLRRM draft manual (February 2001)   
    - GLERM TM 66 ""
    """
    
    #--------------------------------------------------------------------------
    # UNPACK THE INPUT VARIABLES
    #--------------------------------------------------------------------------
    
    # beginning of increment levels
    mhboi = float(BOI_MHU)
    scboi = float(BOI_STC)
    erboi = float(BOI_ERI)
    cgboi = float(BOI_CGI) # to replace niz variable in Matt's code.

    # net total supply for michigan-huron and net basin supply for st.clair and erie
    mhnts = float(NTS_MHU)
    scnbs = float(NBS_STC)
    ernbs = float(NBS_ERI)
    
    # ice and weed retardation values for michigan-huron, st.clair, and erie
    mhice = float(ICW_MHU)
    scice = float(ICW_STC)
    erice = float(ICW_ERI)
        
    # diversion for michigan-huron and erie
    if not DIV_MHU:
        DIV_MHU = 0.0
    if not DIV_ERI:
        DIV_ERI = 0.0
    
    mhdiv = float(DIV_MHU)
    erdiv = float(DIV_ERI)
    
    # consumptive use for michigan-huron, st.clair and erie
    if not CSU_MHU:
        CSU_MHU = 0.0
    if not CSU_STC:
        CSU_STC = 0.0
    if not CSU_ERI:
        CSU_ERI = 0.0
        
    mhcsu = float(CSU_MHU)
    sccsu = float(CSU_STC)
    ercsu = float(CSU_ERI)
    
    # groundwater for michigan-huron, st.clair, and erie
    if not GDW_MHU:
        GDW_MHU = 0.0
    if not GDW_STC:
        GDW_STC = 0.0
    if not GDW_ERI:
        GDW_ERI = 0.0
    
    mhgdw = float(GDW_MHU)
    scgdw = float(GDW_STC)
    ergdw = float(GDW_ERI)       
 
    # seconds per increment (integer)
    dt = deltat

    # constants (float)
    half = 0.5
   
	# outflow coefficients (from midlakes configuration file)
    mhk  = float(MiddleLakesInfo.stagefall_mhk)
    mhym = float(MiddleLakesInfo.stagefall_mhym)
    mha  = float(MiddleLakesInfo.stagefall_mha)
    mhb  = float(MiddleLakesInfo.stagefall_mhb)
    mhwt = float(MiddleLakesInfo.stagefall_mhwt)
    mhc  = float(MiddleLakesInfo.stagefall_mhc)
	
    sck  = float(MiddleLakesInfo.stagefall_sck)
    scym = float(MiddleLakesInfo.stagefall_scym)
    sca  = float(MiddleLakesInfo.stagefall_sca)
    scb  = float(MiddleLakesInfo.stagefall_scb)
    scwt = float(MiddleLakesInfo.stagefall_scwt)
    scc  = float(MiddleLakesInfo.stagefall_scc)
	
    erk  = float(MiddleLakesInfo.stagefall_erk)
    erym = float(MiddleLakesInfo.stagefall_erym)
    era  = float(MiddleLakesInfo.stagefall_era)
    erb  = float(MiddleLakesInfo.stagefall_erb)
    erwt = float(MiddleLakesInfo.stagefall_erwt)
    erc  = float(MiddleLakesInfo.stagefall_erc)
    
    OneHundred = 100
    FlowRound = int(MiddleLakesInfo.roundflo) # -1
    LevelRound = int(MiddleLakesInfo.roundlev) # 4
    
    toleranceMH = float(MiddleLakesInfo.iterativeSCRtolerance)
    toleranceSC = float(MiddleLakesInfo.iterativeDETtolerance)
    toleranceER = float(MiddleLakesInfo.iterativeNIAtolerance)
    
    # Initiate the Converge variable (logical)
    # TODO verify that we should initialize this as True (...False makes more sense to Zoe...)
    Converge = True
    
    # COMPLETED: UNPACKING THE INPUT VARIABLES
    
    #--------------------------------------------------------------------------
    # INITIATE VARIABLES FOR FIRST ITERATION (K=1) "FIRST GUESS OF SOLUTION"
    #--------------------------------------------------------------------------
    
    """ general format of this iterative procedure is to fill in the table:
     
     Iteration |  eoi_guess  |  mlv(eoi_guess) |    q(mlv)  |   eoi_next   
     K-2       |   eoi2ago   |                 |
     K-1       |   eoi1ago   |      imlv       |      q     |   eoi
     K         |   eoi       |
    
    """ 
    
    # define a bogus value (float) for the intial EOI guesses of the Z.
    initial_q_guess = -999.99 # cms [cglrrm -999.99] 4999.99
    initial_m_guess = -999.99 # m
    #initial_dz_guess = 0.50 # 50 cm [initial levles as -999.99 NOT working]
    
    # define the first guess as equal to the boi levels plus some deltaz meters
    mheoi = initial_m_guess # mhboi + initial_dz_guess
    sceoi = initial_m_guess # scboi + initial_dz_guess
    ereoi = initial_m_guess # erboi + initial_dz_guess  
    cgeoi = initial_m_guess # cgboi + initial_dz_guess
    
    # initiate the end of increment level from the previous iteration as float -999.99
    mheoi1ago =  initial_m_guess # mhboi - initial_dz_guess
    sceoi1ago =  initial_m_guess # scboi - initial_dz_guess
    ereoi1ago =  initial_m_guess # erboi - initial_dz_guess   
    
    # initiate the imlv (float) for first outflow calculations
    mhimlv = mhboi
    scimlv = scboi
    erimlv = erboi
    cgimlv = cgboi
    
    # initiate a guess of the mean outflows over the previous increment (outla in cglrrm)
    qmh = initial_q_guess
    qsc = initial_q_guess
    qer = initial_q_guess
    
    # define the maximum number of iterations K (hardcoded limit)
    maxiterations = 9999
    
    print("boi            | %7.2f %7.2f %7.2f" % (mhboi,scboi,erboi))
    print("eoi & flow & e |    mhu     stc     eri |  src    det    nia | e_scr e_det e_nia |")
    
    #-------------------------------------------------------------------------
    # ITERATION CALCULATIONS (1:9999)
    #-------------------------------------------------------------------------
      
    for K in range(1,maxiterations):
        #----------------------------------------------------------------------        
        # define eoi levels as whatever it was for the previous iteration
        #      when K=1, this will be equal to the initalized eoi1ago
        #---------------------------------------------------------------------- 

        # mheoi
        # sceoi
        # ereoi
       
        #----------------------------------------------------------------------  
        # define any eoi2ago100 and eoi1ago and q1ago values
        #----------------------------------------------------------------------  
        
        mheoi2ago100 = int(mheoi1ago*OneHundred)
        sceoi2ago100 = int(sceoi1ago*OneHundred)
        ereoi2ago100 = int(ereoi1ago*OneHundred)
        
        mheoi1ago = mheoi
        sceoi1ago = sceoi
        ereoi1ago = ereoi
        
        qmh1ago = qmh
        qsc1ago = qsc
        qer1ago = qer        
       
        # --------------------------------------------------------------------
        # calculate what flows for the increment would be if eoi is above
        # --------------------------------------------------------------------
        
        outmh = Outflow(mhimlv,scimlv,mhk,mhym,mha,mhb,mhwt,mhc)
        outsc = Outflow(scimlv,erimlv,sck,scym,sca,scb,scwt,scc)
        outer = Outflow(erimlv,cgimlv,erk,erym,era,erb,erwt,erc)
        
        qmh = outmh - mhice
        qsc = outsc - scice
        qer = outer - erice    
        
        qmh = round(qmh, FlowRound)
        qsc = round(qsc, FlowRound)
        qer = round(qer, FlowRound)

        # --------------------------------------------------------------------
        # calculate the eoi value based on the eoi guess for this iteration
        # --------------------------------------------------------------------
        # NOTE: 
        # - if the eoi guess above was the true solution, then the new eoi 
        #        value calculated below would be equal to the guess
        #        (redefined as eoi1ago prior to the new eoi calculation) 
        # - if the eoi guess above was not the true solution, then the new eoi
        #        value calculated below would not be equal to the guess
        #        (redefined as eoi1ago prior to the new eoi calculation) 
        
        Amh = GetArea('mh',mhimlv)
        Asc = GetArea('sc',scimlv)
        Aer = GetArea('er',erimlv)
        
        mheoi = mhboi + ((mhnts       + mhgdw) - (mhdiv + mhcsu + qmh)) * dt/Amh
        sceoi = scboi + ((scnbs + qmh + scgdw) - (        sccsu + qsc)) * dt/Asc
        ereoi = erboi + ((ernbs + qsc + ergdw) - (erdiv + ercsu + qer)) * dt/Aer
        
        #mheoi = mhboi + (mhnts - mhdiv       - qmh) * dt/Amh
        #sceoi = scboi + (scnbs         + qmh - qsc) * dt/Asc
        #ereoi = erboi + (ernbs - erdiv + qsc - qer) * dt/Aer
        
        mheoi = round(mheoi,LevelRound)
        sceoi = round(sceoi,LevelRound)
        ereoi = round(ereoi,LevelRound)
        
        #----------------------------------------------------------------------
        # calculate what the mean value for the increment would be if eoi is above
        #----------------------------------------------------------------------
        # NOTE - original code has mean levels calculated after outflows and eois...?
        mhimlv = half * (mhboi + mheoi)
        scimlv = half * (scboi + sceoi)
        erimlv = half * (erboi + ereoi)
        cgimlv = half * (cgboi + cgeoi)                        
        
        #----------------------------------------------------------------------
        # calculate the error from the flow calculated for this iteration
        # and the previous iteration
        #----------------------------------------------------------------------
        
        diffQmh = abs(qmh - qmh1ago)
        diffQsc = abs(qsc - qsc1ago)
        diffQer = abs(qer - qer1ago)
        
        # ---------------------------------------------------------------------
        # print for logging[zoe working on this part...to "see" iterations]
        # ---------------------------------------------------------------------
        print("iteration K %3d %7.2f %7.2f %7.2f %6.0f %6.0f %6.0f %7d %5.0f %5.0f"
              % (K, mheoi1ago, sceoi1ago, ereoi1ago, qmh, qsc, qer, diffQmh, diffQsc, diffQer))
                    
        #----------------------------------------------------------------------
        """ SCENARIO 1: CONVERGE ... flow convergence criteria met ... """
        #----------------------------------------------------------------------
        # What does this mean? If the |Q_K - Q_Kprior| <= QTolerance for all
        # channels then the solution has converged towards the true solution 
        # and we can stop iterating.  If just one connecting channel does not 
        # meet the criteria, then this tolerance is not met.
        
        if diffQmh <= toleranceMH and diffQsc <= toleranceSC and diffQer <= toleranceER:

            Converge = True
           
        else:
            # however, before going through the rest of the iterations
            # let's check a partiuclar scenario that might occur which is
            # outlined below. TODO therefore no break here? keep going...
            Converge = False
 
        #----------------------------------------------------------------------
        """ SCENARIO 2: CONVERGE ... oscillating levels criteria met ... """
        #----------------------------------------------------------------------
        # What does this mean? Well, maybe the flow criteria is not being met.
        # However, consider that maybe for each iteration we are bouncing
        # around the true solution, just mearly outside the tolerance defined.
        
        # Example: K = 37  ,  eoi =  175.87    diffQ = 11    (tolerance = 10)
        #          K = 38  ,  eoi =  175.85    diffQ = 12    (tolerance = 10)
        #          K = 39  ,  eoi =  175.87    diffQ = 11    (tolerance = 10)
        #          K = 40  ,  eoi =  175.85    diffQ = 12    (tolerance = 10)
        #          K = 41  ,  eoi =  175.87    diffQ = 11    (tolerance = 10)        
                 
        if LevelRound != -99:
            # NOTE: round(173.32412,-99) = 0.0
            
            # OneHundred is 100.
            # int(eoi*100) ~ int(175.87333 * 100) ~ 17587

            # TODO ZM review original code make sure nothing missing or wrong!
            
            if mheoi2ago100 == int(mheoi*OneHundred) and sceoi2ago100 == int(sceoi*OneHundred) and ereoi2ago100 == int(ereoi*OneHundred):
                
                # mheoi2ago = int(mheoi*100) ~ 17687
                # sceoi2ago = int(sceoi*100) ~ 17554 
                # ereoi2ago = int(ereio*100) ~ 17467    
     
               Converge = True
               
               # TODO make sure convergenec makes sense - if first convergence
               # criteria is met then does the second always get med? else - break?
                
            else:
                                             
                Converge = False
                     
        #----------------------------------------------------------------------
        # IF THE SOLUTION DOES NOT CONVERGE UNDER EITHER CRITERIA
        #----------------------------------------------------------------------
        
        # the following variables are used in the next iteration
        
        # eoi: this guess of eoi is used in the next iteration
        # eoi1ago: this guess is used to tabulate eoi2ago100 in next iteration
        # out: this approximation of out will be used to define out1ago in next iteration
        
        # go to next iteration!

        #----------------------------------------------------------------------
        # IF THE SOLUTION CONVERGES UNDER AT LEAST ONE CRITERIA
        #----------------------------------------------------------------------
        
        if Converge == True:
            # stop iterating!
            # accept the current eoi levels as the solution
            break

    #--------------------------------------------------------------------------
    # SUMMARY OF END OF INCREMENT
    #--------------------------------------------------------------------------
   
    if Converge == True:
        # how many iterations did we go through for this increment?
        # print(K)
        
        return round(mheoi,6),round(sceoi,6),round(ereoi,6)
# IMPORTANT: THIS CODE IS INCOMPLETE AND UNOFFICIAL AND WORK IN PROGRESS!!!!!
#  - sharing with group just so everyone can see how things are evolving
# import sys
import datetime

#------------------------------------------------------------------------------
#                          ABOUT: midlakes_util.py
#------------------------------------------------------------------------------
"""utilities requires for running midlakes
CLASS: ml_config_info()
       this will turn a config file into an object of class ml_config_info
       
FUNCTIONS: set of functions required for created the ml_config_info class
      - line_for_parsing(line)
      - parse1(line)
      - FindLineString(ConfigHeader,filename)
      - FindLineDate(ConfigDate,filename)
      - FindLineFloat(ConfigFloat,filename)
      - FindLineInt(ConfigFloat,filename)
      - FindLineOutflow(ConfigOutflowParam,filename)

      - DaysinPeriod(PeriodString)
      - SecondsinPeriod(PeriodString)
      - SecondsinIncrement(SecondsInPeriod,NumIncrements)
      - MiddleLakesSolution(SolutionString)
      - FormatRoutingKey(PeriodString)
      
FUNCTIONS: for midlakes program
      - GetArea(LAKENUM, Z)
      - Outflow(ZU, ZD, K, YM, A, B, WT, ZC)
      
      - GetLevelAtUSSlip(StMarysFlow,MHLevel,USSlipEqArg,MonthIndex,UseDefaultAdj,objParam):
      - USSlipParam (OBJECT)
"""

#----------------------------------------------------------------------------
# Class object full of inputs for midlakes routing
#----------------------------------------------------------------------------

class ml_config_info():
    def __init__(self,filename):
        self.filename = filename
        """
        Object full of the middle lakes configuration file information
		
		- configuration file
		- title
        """
        self.model = FindLineString('master model',filename)
        self.titleI = FindLineString('model title i',filename)
        self.titleII = FindLineString('model title ii',filename)
        self.sdate = FindLineDate('startdate',filename)
        self.edate = FindLineDate('enddate',filename)
        self.startmhu = FindLineFloat('starting mhu level',filename)
        self.startstc = FindLineFloat('starting stc level',filename)
        self.starteri = FindLineFloat('starting eri level',filename)
        self.outdir = FindLineString('output directory',filename)
        self.mlsolution = FindLineString('solution method',filename) # m --> 1, i --> 2
        self.units = FindLineString('units',filename) # English or Metric
        self.mlrouting = FindLineString('routing time step',filename)
        self.mlincrements = FindLineInt('increments',filename)
        self.mlverbosity = FindLineInt('verbosity',filename)
        self.iterativeSCRtolerance = FindLineFloat('tolerance scr',filename)
        self.iterativeDETtolerance = FindLineFloat('tolerance det',filename)
        self.iterativeNIAtolerance = FindLineFloat('tolerance nia',filename)
        self.timedelta_one_day = datetime.timedelta(1)
        
        # Routing Information
        # secinroutingperiod
        # secinincrement
        # daysinroutingperiod
        # daysinroutingperiodcount
        # keyintvl
        self.keyintvl = FormatRoutingKey(self.mlrouting)
		
        # Thresholds for Checks (conservation of mass)
        self.roundlev = FindLineInt('rounding levels',filename)
        self.roundflo = FindLineInt('rounding flows',filename)
        
        # Middle Lakes Solution
        self.solutionindicator = MiddleLakesSolution(self.mlsolution)
        
        # Outflow Coefficients
        self.stagefall_mhk = FindLineOutflow('outflow parameters scr',filename)[0]
        self.stagefall_mhym = FindLineOutflow('outflow parameters scr',filename)[1]
        self.stagefall_mha = FindLineOutflow('outflow parameters scr',filename)[2]
        self.stagefall_mhb = FindLineOutflow('outflow parameters scr',filename)[3]
        self.stagefall_mhwt = FindLineOutflow('outflow parameters scr',filename)[4]
        self.stagefall_mhc = FindLineOutflow('outflow parameters scr',filename)[5]

        self.stagefall_sck = FindLineOutflow('outflow parameters det',filename)[0]
        self.stagefall_scym = FindLineOutflow('outflow parameters det',filename)[1]
        self.stagefall_sca = FindLineOutflow('outflow parameters det',filename)[2]
        self.stagefall_scb = FindLineOutflow('outflow parameters det',filename)[3]
        self.stagefall_scwt = FindLineOutflow('outflow parameters det',filename)[4]
        self.stagefall_scc = FindLineOutflow('outflow parameters det',filename)[5]

        self.stagefall_erk = FindLineOutflow('outflow parameters nia',filename)[0]
        self.stagefall_erym = FindLineOutflow('outflow parameters nia',filename)[1]
        self.stagefall_era = FindLineOutflow('outflow parameters nia',filename)[2]
        self.stagefall_erb = FindLineOutflow('outflow parameters nia',filename)[3]
        self.stagefall_erwt = FindLineOutflow('outflow parameters nia',filename)[4]
        self.stagefall_erc = FindLineOutflow('outflow parameters nia',filename)[5]

        # Print the config file

#-----------------------------------------------------------------------------
# Functions for parsing [Written by Tim Hunter]
#------------------------------------------------------------------------------
#  Given a line of text (i.e. a string)
#  1. strip off any comments (everything at/after the first # character)
#  2. convert to all lowercase
#
def line_for_parsing(line):
    """line_for_parsing(line)
    ABOUT: Parse the line (ignore #'s)"""
    i = line.find('#')
    s = line
    if (i >= 0):
       s = line[0:i-1].rstrip()
    return s.lower()


#---------------------------------------------------------------------------------
#  Given a line of text (i.e. a string)
#  Look for the first colon.  If you find one return a tuple with
#  everything to the left of the colon and everything to the right
#  of the colon.
#  If no colon is found, return None.
#  e.g.
#     'startdate: 2001,01,01'          -> ('startdate', ' 2001,01,01')
#     'starttime: 2001-01-01 18:00:00' -> ('starttime', ' 2001-01-01 18:00:00')
#     'startdate= 2001,01,01'          -> None
#
def parse1(line):
    i=line.find(':')
    if (i > 0):
        a = line[0:i]
        b = line[i+1:]
        return a,b
    
#------------------------------------------------------------------------------
# Functions to find and read in the header information [Written by Zoe Miller]
#------------------------------------------------------------------------------
def FindLineString(ConfigHeader,filename):
    """FindLineString(ConfigHeader,filename)
    ABOUT: Fine the string info for line starting with ConfigHeader"""
    Header = ConfigHeader
    HeaderInfo = None 
    
    with open(filename,"r") as f:
        for line in f:
            s1 = line_for_parsing(line)
            if s1.find('#') < 0:
                p=parse1(s1)
                if (p):
                    if p[0].strip() == Header:
                        s = p[1].strip()
                        HeaderInfo = s
    return HeaderInfo

def FindLineDate(ConfigDate,filename):
    Header = ConfigDate
    adate = None 
    
    with open(filename,"r") as f:
        for line in f:
            s1 = line_for_parsing(line)
            if s1.find('#') < 0:
                p=parse1(s1)
                if (p):
                    if p[0].strip() == Header:
                        y,m,d = p[1].split(',')
                        adate = datetime.date(int(y), int(m), int(d))
    return adate

def FindLineFloat(ConfigFloat,filename):
    Header = ConfigFloat
    ValueFloat = None 
    
    with open(filename,"r") as f:
        for line in f:
            s1 = line_for_parsing(line)
            if s1.find('#') < 0:
                p=parse1(s1)
                if (p):
                    if p[0].strip() == Header:
                        s = p[1].split()
                        ValueFloat = float(s[0])
    return ValueFloat

def FindLineInt(ConfigInteger,filename):
    Header = ConfigInteger
    ValueInt = None 
    
    with open(filename,"r") as f:
        for line in f:
            s1 = line_for_parsing(line)
            if s1.find('#') < 0:
                p=parse1(s1)
                if (p):
                    if p[0].strip() == Header:
                        s = p[1].split()
                        ValueInt = int(s[0])
    return ValueInt


def FindLineOutflow(ConfigOutflowParam,filename):
    """FindLineOutflow(ConfigOutflowParam,filename)
    INPUTS: Header of format "outflow parameters riv
    OUTPUT: ARRAY/VECTOR of [k,ym,a,b,wt,c]"""
    Header = ConfigOutflowParam
    k = None
    ym = None 
    a = None 
    b = None 
    wt = None
    c = None 
    
    with open(filename,"r") as f:
        for line in f:
            s1 = line_for_parsing(line)
            if s1.find('#') < 0:
                p=parse1(s1)
                if (p):
                    if p[0].strip() == Header:
                        s = p[1].split()
                        k = float(s[0])
                        ym =  float(s[1])
                        a =  float(s[2])
                        b =  float(s[3])
                        wt =  float(s[4])
                        c =  float(s[5])
    return k, ym, a, b, wt, c
        
# ----------------------------------------------------------------------------
# Functions to transform middlakes configuration information to program criteria         
# ----------------------------------------------------------------------------

def MiddleLakesSolution(SolutionString):
    SecondsPerIncrement = SolutionString.lower()
    FirstChar = SecondsPerIncrement[0]
    if FirstChar == 'i':
        # 'Iterative','iterative solution' or anything that starts with I or i.
        Solution = 1
    elif FirstChar == 'm':
        # 'Matrix','matrix sol' or anything that starts with M or m.
        Solution = 2
    elif FirstChar == 'n':
        Solution = 3
    return Solution

def FormatRoutingKey(PeriodString):
    PeriodLower = PeriodString.lower()
    FirstChar = PeriodLower[0]
    if FirstChar == 'd':
        FormatRoutingKeyName = 'dy'
        return FormatRoutingKeyName
    elif FirstChar == 'w':
        FormatRoutingKeyName = 'wk'
        return FormatRoutingKeyName    
    elif FirstChar == 'q':
        FormatRoutingKeyName = 'qm'
        return FormatRoutingKeyName     
    elif FirstChar == 'm':
        FormatRoutingKeyName = 'mn'
        return FormatRoutingKeyName

#----------------------------------------------------------------------------
#                             GetArea function
#----------------------------------------------------------------------------

def GetArea(LAKENUM = None, Z = None):
    """
    - LAKENUM {1,2,5,6,7} for {Sup,MiH,StC,Eri,Ont} respectively
    - Z (WLELEVATION - not used here and wasn't "really" used in original code)
    - GETAREA returns area in square meters
    
    Original CGLRRM script had a GETAREA function that produced an area that
    was dependent on lake level.  In reality this is true -- the lake area does
    change with respect to level, however in previous operational set-up of the 
    CGLLRM the area was always assumed to be constant, irregardless of level. 
    This is standard procedure/assumption for this scope of hydrologic modeling. 
    """
    
    if not LAKENUM:
        raise Exception('Lake number not specified')
    if not Z:
        raise Exception('No lake leve specified...Required although not used')
        
    if LAKENUM=='sp':     # Lake Superior        ~  82100.00 # km2
        AREA=82100000000.00 # sq m 
    elif LAKENUM=='mh':   # Lakes Michigan-Huron ~ 117400.00 # km2
        AREA=117400000000.00 # sq m
    elif LAKENUM=='sc':   # Lake St. Clair       ~   1114.00 # km2
        AREA=1114000000.00 # sq m          
    elif LAKENUM=='er':   # Lake Erie            ~  25700.00 # km2
        AREA=25700000000.00 # sq m         
    elif LAKENUM=='on':   # Lake Ontario         ~  18960.00 # km2
        AREA=18960000000.00 # sq m  
    else:
        AREA=0.0
    return float(AREA)

#----------------------------------------------------------------------------
#                           Outflow function
#----------------------------------------------------------------------------

def Outflow(ZU = None, ZD = None , K = None, YM = None, A = None, B = None, WT = None, ZC = None):
    """ FUNCTION USED IN MIDLAKES:
        
    INPUTS REQUIRED:
	
    zu = upstream lake elevation (m)
    zd = downstream lake elevation (m)
    k = outflow equation coefficient
    ym = outlet bottom elevation
    a = depth exponent
    b = fall exponent
    wt = relative weight between upstream and downstream lake levels
    zc = fall constant 
   
    OUTPUTS: OUTFLOW (cms)
    
    ABOUT:
        WTDEPTH = WT*ZU + (ONE-WT)*ZD - YM
        FALL = ZU-ZD+ZC
        OUTFLOW = K*(WTDEPTH)**A*(FALL)**B
    
        NOTE: When b is equal to zero, the resulting equation is single-gage or
        stage discharge equation. Flow will be independent of downstream level
        if wt is equal to one.
        
    CODE HISTORY:
        
    This function was translated from FORTRAN into python by Zoe Miller (USACE)
    in December 2018.  A slight modification to only run through the feasability
    checks when b is not equal to zero, instead of always.

    Previous code was written by Matt McPherson (previously USACE, now HEC). He
    added the ZC variable to accomodate the metric forms of the equations. [ZM:
    when was this term added...any sources?]
    
    Matt references code originally written by Anne H. Clites, NOAA/GLERL, 
    September, 1995.  The constants k, ym, a, and b are dictated by the 
    stage-discharge relationship for each connected channel. """
    
    if not isinstance(ZU,float):
        raise Exception('No upstream elevation (zu) specified')
    if not isinstance(ZD,float):
        raise Exception('No downstream elevation (zd) specified')
    if not isinstance(K,float):
        raise Exception('No outflow coefficient (k) specified')
    if not isinstance(YM,float):
        raise Exception('No outlet bottom elevation (ym) specified')
    if not isinstance(A,float):
        raise Exception('No depth exponent (a) specified')
    if not isinstance(B,float):
         raise Exception('No fall exponent (b) specified')
    if not isinstance(WT,float):
        raise Exception('No weigt (wt) specified')
    if not isinstance(ZC,float):
        raise Exception('No fall constant (zc) specified') 
    
    # equation inputs
    # ZU, ZD, K, YM, A, B, WT, ZC
    
    # constant
    one = 1.0
    
    #-----------------------------------------------------------------------
    # CALCULATE THE FLOW
    #-----------------------------------------------------------------------
    
    # Only run through check A if two-gage equation (b not equal to zero).
    if B != 0:
        if ZD >= (ZU+ZC):
            # FALL TERM VALIDITY CHECK
            OUTFLOW = 'NA'
            print('downstream level is higher than upstream! : ZD',ZD,'>=','ZU+ZC',ZU+ZC)
            return OUTFLOW
    
    # Run through check B for all scenarios    
    if YM > ZU:
        OUTFLOW = 'NA'
        print('outlet bottom elevation is higher than upstream!: YM',YM,'>','ZU',ZU)
        return OUTFLOW
      
    """
    WTDEPTH = WT*ZU + (one-WT)*ZD - YM
    FALL = ZU-ZD+ZC
    OUTFLOW = K*(WTDEPTH**A)*(FALL**B)
    """
    
    WTDEPTH = WT*ZU + (one-WT)*ZD - YM
    FALL = ZU-ZD+ZC
    OUTFLOW = K*(WTDEPTH**A)*(FALL**B)

    return OUTFLOW

#b1 = 1.0
#ze = 0.0

# TEST THE OUTFLOW FUNCTION

#d = Outflow(174.2, 0, 669.4, 170.043, 1.5, b1,1.0,ze) 


#------------------------------------------------------------------------------
# Write 
#------------------------------------------------------------------------------


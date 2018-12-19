# IMPORTANT: THIS CODE IS INCOMPLETE AND UNOFFICIAL AND WORK IN PROGRESS!!!!!
#  - sharing with group just so everyone can see how things are evolving
#------------------------------------------------------------------------------
#                          Middle Lake Solution: Matrix
#------------------------------------------------------------------------------
# Coded in python by Zoe Miller zoe.a.miller@usace.army.mil
# Last Edit: 11/8/2018
# 
# BRIEF OVERVIEW: 
# This is a python script contains one function: SolveMatrix
# This is adapted from the subroutine "CoeffMatrix" written by Matt McPherson
# for the CGLRRM in the midlakes.for script in FORTRAN. His code was an 
# adaptation from midlakes.for code by Anne H. Clites, NOAA/GLERL, Nov 1995.
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
#            How to Solve the Problem: Matrix Solution
# ---------------------------------------------------------------------------
# ABOUT MATRIX SOLUTION:  

# 1) Recall that NBS, D, R & QO's are rates over the time step (dt)
# 2) Represent the connecting channel flot terms as projecting half way into 
#    the time step.  Think: y = mt + b = b + mt.  Recall m is the slope.
#    Now think recall that derivative is best estimate of slope evaluated at t
#              f(x(t+dt)) = f(x) + f'(x(t))*dt
#    Now note that Q is a function of two variables. And project halfway into dt
#          Q(z1(t+dt/2),z2(t+dt/2)) = Q(z1(t),z2(t)) + (dt/2)(Q'(z1(t),z2(t)) 
# 3) Substitution: Replace the QO's = QOt + (1/2)(dQ/dt)dt [**KEY TRICK**]
#        - See Quinn (1978) and Midlakes NOAA GLERL TM
# 4) Calculate and evaluate the partials at the beginning of the time step
#     and calculate the flows at the beginning of the time step for QO1.
# 5) Substutite in all the terms calculated. Observe the 3 equations become
#    a system of which the only unknowns are dzs.  
# 6) Rearrange the knowns and unknowns into a format of Ax=b.
# 7) Solve the linear system: Use Cramer's rule to solve the system of equations
#        - http://mathworld.wolfram.com/CramersRule.html 
# 8) The method above solves for dzs. Final step is to calculate z_{n}
#        -z_{n+1} = z{n+1} + dz
#  
# -----------------------------------------------------------------------------
#                 Import python modules/packages or functions
#------------------------------------------------------------------------------

from middle.routemiddle_util import Outflow, GetArea 

#----------------------------------------------------------------------------
#                         Function: SolveMatrix
#----------------------------------------------------------------------------

def SolveMatrix(BOI_MHU=None,BOI_STC=None,BOI_ERI=None,BOI_CGI=None,
                NTS_MHU=None,NBS_STC=None,NBS_ERI=None,
                ICW_MHU=None,ICW_STC=None,ICW_ERI=None,
                DIV_MHU=None             ,DIV_ERI=None,
                CSU_MHU=None,CSU_STC=None,CSU_ERI=None,
                GDW_MHU=None,GDW_STC=None,GDW_ERI=None,
                deltat = None,
                MiddleLakesInfo=None):
    
    """
    Middle Lakes Solution: Matrix Solution
    
    ABOUT: 
    Function solves the "hydrologic equations" for one increment interval using
    the matrix solution.
    
    INPUTS:
    12 float variables: 
    boi: mhu, stc, eri, cgi
    nts: mhu
    nbs: stc, eri
    icw: mhu, stc, eri
    div: mhu, eri
    & 1 variable object of class/type "ml_config_inf". see middle_util.py for more details.
    MiddleLakesInfo
    """
    #-----------------------------------------------------------------------
    # UNPACK VARIABLES NEEDED FOR SOLUTION CALCULATIONS
    #-----------------------------------------------------------------------
    
    # "Unpack" the 12 floats
    mhboi = BOI_MHU
    scboi = BOI_STC
    erboi = BOI_ERI
    cgboi = BOI_CGI
    # cgboi =  171.03 ZM - where in the world did I get this number from?
    
    # net basin/total supplies during the time step
    mhnts = NTS_MHU
    scnbs = NBS_STC
    ernbs = NBS_ERI
    
    # ice and weed retardation values during the time step
    mhice = ICW_MHU
    scice = ICW_STC
    erice = ICW_ERI
    
    # diversions during the time step
    mhdiv = DIV_MHU
    erdiv = DIV_ERI

    # define constants
    zero = 0.0
    half = 0.5
    one = 1.0

    # UNPACK the seconds per increment (dt)
    dt = deltat
    
    # UNPACK the level rounding
    roundlev = int(MiddleLakesInfo.roundlev) 
    
    # UNPACK stage-fall outflow coefficients
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

    #-----------------------------------------------------------------------
    # BEGIN CALCULATIONS
    #-----------------------------------------------------------------------

    
    #-----------------------------------------------------------------------
    # DETERMINE TWO-GAGE STAGE-FALL OUTFLOW FOR PARTIAL DERIVATIVES FOR COEFFICIENTS
    #-----------------------------------------------------------------------
    # note all shall be evaluated at the beginning of the interval
    
    # outmh :(ice/weed free) outflow of michigan-huron (st.clair river)
    # outmh1:(ice/weed free) outflow of scr evaluated when b ~ b-1
    # outmh2:(ice/weed free) outflow of scr evaluated when a ~ a-1
    
    # outsc :(ice/weed free) outflow of st. clair (detroit river)
    # outsc1:(ice/weed free) outflow of det evaluated when b ~ b-1
    # outsc2:(ice/weed free) outflow of det evaluated when a ~ a-1
    
    # outer :(ice/weed free) outflow of erie (nia river)
    # outer1:(ice/weed free) outflow of nia evaluated when b ~ b-1
    # outer2:(ice/weed free) outflow of nia evaluated when a ~ a-1
    
    # OUTFLOW FCN   | Zup | Zdown|  K |  ym |     a    |     b    |  wt |  c |
    outmh  = Outflow(mhboi, scboi, mhk, mhym, mha      , mhb      , mhwt, mhc)
    outmh1 = Outflow(mhboi, scboi, mhk, mhym, mha      , mhb - one, mhwt, mhc)
    outmh2 = Outflow(mhboi, scboi, mhk, mhym, mha - one, mhb      , mhwt, mhc)

    outsc  = Outflow(scboi, erboi, sck, scym, sca      , scb      , scwt, scc)
    outsc1 = Outflow(scboi, erboi, sck, scym, sca      , scb - one, scwt, scc)
    outsc2 = Outflow(scboi, erboi, sck, scym, sca - one, scb      , scwt, scc)

    outer  = Outflow(erboi, cgboi, erk, erym, era      , erb      , erwt, erc)
    outer1 = Outflow(erboi, cgboi, erk, erym, era      , erb - one, erwt, erc)
    outer2 = Outflow(erboi, cgboi, erk, erym, era - one, erb      , erwt, erc)
    # NIAGARA NOTE: when niaraga b (erb) is zero, both outer and outer2 become
    # single-gage which means that cgboi does not eactually impact the outflow
    # (but is required as an input, irregardles). 

    #-----------------------------------------------------------------------
    # DETERMINE THE COMPONENTS OF THE PARTIALS
    #-----------------------------------------------------------------------
    # partial derivate overview
    #        dQ/dt = (dQ/dz1)(dz1/dt) + (dQ/dz2)(dz2/dt)
    #
    # half way into time step amd evaluated with boi
    #     (dt/2)*dQ/dt = (dt/2)(dQ/dz1)(dz1/dt) + (dQ/dz2)(dz2/dt)
    #     (dt/2)*dQ/dt = (1/2)(dQ/dz1)dz1 + (1/2)(dQ/dz2)dz2
    #
    # Evaluate partials at boi and cut in half.
    # (1/2)(dQ/dz1)
    # (1/2)(dQ/dz2)
    pscrdzmh = half*( mhb*outmh1 + mha*     mhwt *outmh2)
    pscrdzsc = half*(-mhb*outmh1 + mha*(one-mhwt)*outmh2)
    
    pdetdzsc = half*( scb*outsc1 + sca*     scwt *outsc2)
    pdetdzer = half*(-scb*outsc1 + sca*(one-scwt)*outsc2)
    
    pniadzer = half*( erb*outer1 + era*     erwt *outer2)
    pniadzcg = half*(-erb*outer1 + era*(one-erwt)*outer2)
    # NIAGARA NOTE: As explained above, when b=0, outer and outer2 are zero.
    # Therefore the dz4 (see below) multiplier of pniadzcg is actually zero.
    #      (-erb*outer1 + era*(one-erwt)*outer2)  = 0
    #      (   0*outer1 + era*(one-erwt)*0)       = 0   
    #                                    pniadzcg = 0
    
    #-----------------------------------------------------------------------
    # CALCULATE THE COEFFICIENTS OF DZ TERMS OR TRUELY MATRIX 'A' of Ax=b
    #-----------------------------------------------------------------------    
 
    # NOTE: clgrrm had GetArea lakes id by a # 2,5,6 for "mh","sc","er" 
    # Define the Areas
    Amh = GetArea("mh",mhboi)
    Asc = GetArea("sc",scboi)
    Aer = GetArea("er",erboi)
    
    c11=  pscrdzmh + Amh/dt
    c12=  pscrdzsc
    c13=  zero
    
    c21= -pscrdzmh
    c22=  pdetdzsc - pscrdzsc + Asc/dt
    c23=  pdetdzer
    
    c31=  zero
    c32= -pdetdzsc
    c33=  pniadzer - pdetdzer + Aer/dt
    
    #-----------------------------------------------------------------------
    # DETERMINE WATERBALANCE OR RIGHT HAND SIDE OF MATRIX 'b' of Ax=b
    #-----------------------------------------------------------------------
    
    # Sum the known water balance (wb) components: NBS, diversions, consumptive use,
    # inflows, outflows and ice and weed retardations.  These wb values will 
    # become the RHS('b') values of a matrix problem of the form Ax=b.
    
    # cgip dz
    dz4 = 0.0
    
    # water balance terms
    mhwb = mhnts - mhdiv - (outmh - mhice)
    scwb = scnbs + (outmh - mhice) - (outsc - scice)
    erwb = ernbs - erdiv + (outsc - scice) - (outer - erice) - pniadzcg * dz4

    # What in the world is dz4? Technically it is in equation 9 from TM-109
    # Note that this term is multiplied by a term (pniadzcg) that is zero. 
    # Technically this term arises from taking the partial derivative of the 
    # niagara outflow, when it is a two-gage equation. For the sake of keeping 
    # the equations acurate, although it has no influence on the results, 
    # it is included in the code in the water balance total.
    
    #-----------------------------------------------------------------------
    # REVIEW OF MATRIX PROBLEM: MATRIX OF FORM Ax=b
    #-----------------------------------------------------------------------
    # Status of Problem:
    # 
    # c11*dZ1 + c12*dZ2 + c13*dZ3 = C1 = mhwb
    # c21*dZ1 + c22*dZ2 + c23*dZ3 = C2 = scwb
    # c31*dZ1 + c32*dZ2 + c33*dZ3 = C3 = erwb 
    #
    # |c11 + c12 + c13| |dZ1|   |C1|   |mhwb|
    # |c21 + c22 + c23| |dZ2| = {C2} = {scwb}
    # |c31 + c32 + c33| |dZ3|   |C3|   |erwb|
    
    # The problem has boiled down to linear system that is solvable explicility
    # A system of the form Ax=b can be solved using Cramer's Rule.
    # 
    # The coefficients of cij make up the matrix C (think: "A")
    # The water balence terms make up the vector WB (think: "b")
    # The dz terms make up the unknown vector dz (think: "x")

    #-----------------------------------------------------------------------
    # CALCULATE DETERMINATES TO SOLVE THE MATRIX
    #-----------------------------------------------------------------------
    """    
    This matrix problem is solved by using determinats and Cramer's Rule:
        
        Wikipedia: https://en.wikipedia.org/wiki/Determinant
        Wikipedia: https: https://en.wikipedia.org/wiki/Cramer%27s_rule
        Wolfram/MathWorld: http://mathworld.wolfram.com/CramersRule.html

    Cramer's Rule let's us know that the solution (dZ1,dZ2,dZ3) can be solved by:
    
          |mhwb  c12  c13|         |c11  mhwb  c13|         |c11 + c12  mhwb|
          |scwb  c22  c23|         |c21  scwb  c23|         |c21 + c22  scwb|
          |erwb  c32  c33|         |c31  erwb  c33|         |c31 + c32  erwb|
    dZ1 = ------------------ , dZ1 = ------------------ , dZ1 = ------------------     
          |c11  c12  c13|          |c11  c12  c13|          |c11  c12  c13|
          |c21  c22  c23|          |c21  c22  c23|          |c21  c22  c23|
          |c31  c32  c33|          |c31  c32  c33|          |c31  c32  c33|
    
    And each numerator and denominator can be calcuated using determinats of the
    smaller matrices.        
    """
    
    #  cdet is the determinant of the C matrix         [denominator for all] 
    cdet = c11*((c22*c33)-(c23*c32)) - \
        c12*((c21*c33)-(c23*c31)) + \
        c13*((c21*c32)-(c22*c31))
    
    #  z1det is the determinant for Michigan-Huron.      [numerator for dz1] 
    z1det = mhwb*((c22*c33)-(c23*c32)) - \
        c12*((scwb*c33)-(c23*erwb)) + \
        c13*((scwb*c32)-(c22*erwb))
    
    #  z2det is the determinant for St. Clair.           [numerator for dz2] 
    z2det = c11*((scwb*c33)-(c23*erwb)) - \
        mhwb*((c21*c33)-(c23*c31)) + \
        c13*((c21*erwb)-(scwb*c31))
    
    #  z3det is the determinant for Erie.                [numerator for dz3] 
    z3det = c11*((c22*erwb)-(scwb*c32)) - \
        c12*((c21*erwb)-(scwb*c31)) + \
        mhwb*((c21*c32)-(c22*c31))
    
    # Solved System: DeltaZs
    dz1 = z1det/cdet
    dz2 = z2det/cdet
    dz3 = z3det/cdet

    #-----------------------------------------------------------------------
    # CALCULATE NEW EOI VALUES USING THE CALCULATED DELTA Z's: EOI = BOI + DZ
    #-----------------------------------------------------------------------
    
    mheoi = mhboi + dz1
    sceoi = scboi + dz2
    ereoi = scboi + dz3
     
    #-----------------------------------------------------------------------
    # OUTPUT SHALL BE THE EOI LEVELS
    #-----------------------------------------------------------------------
    return round(mheoi,roundlev), round(sceoi,roundlev), round(ereoi,roundlev)


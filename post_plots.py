
#----------------------------------------------------------------
#  This is a plot to do the scripting for midlakes after it runs CGLRRM
#-----------------------------------------------------------------
# Start of a post processing script for plotting the results
# length of the datetime 

import datetime as dt

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import numpy as np

import handler.databank as db
# ------------ import created modules: glrrm tools? ----------------
# functions used for between superior/middle/ontario? or would this be from handler?

def daterangeBOP(VaultData,ml_info):
    # ------- info from ml_info object------
    #---------------------------------------
    key_intvl = ml_info.keyintvl
    bop = ml_info.sdate
    
    tfirst = ml_info.sdate
    tlast = ml_info.edate
    
    datlen = len(VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='smr', first=tfirst , last=tlast).dataVals)
    # initialize bop list (empty)
    list_bop = []
    # -------------------------------------
    for r in range(0,datlen):
        if key_intvl[0].lower()=='d':
            list_bop.append(bop)
            eop = bop
            lenperiodtd = dt.timedelta(days=1)
            lenperiodc = lenperiodtd.days
            bop = eop + dt.timedelta(days=1)
        if key_intvl[0].lower()=='w':
            list_bop.append(bop)
            eop = bop + dt.timedelta(days=6)
            lenperiodtd = dt.timedelta(days=7)
            lenperiodc = lenperiodtd.days
            bop = eop + dt.timedelta(days=1)
        if key_intvl[0].lower()=='q':
            list_bop.append(bop)
            qtr = db.util.qtr_month_from_date(bop)
            lenperiodc = db.util.days_in_qtr_mon(bop.year,bop.month,qtr)
            lenperiodtd = dt.timedelta(days=lenperiodc)
            eopd = db.util.getQtrMonthStartEnd(bop.year,bop.month,qtr)[1]
            eop = dt.datetime(bop.year,bop.month,eopd)
            lenperiodc = lenperiodtd.days
            bop = eop + dt.timedelta(days=1)
        if key_intvl[0].lower()=='m':
            list_bop.append(bop)
            lenperiodc = db.util.days_in_month(bop.year,bop.month)
            lenperiodtd = dt.timedelta(days=lenperiodc)
            
            eop = db.util.last_day_of_month(bop)
            bop = eop + dt.timedelta(days=1)

    # Create some arrays for the        
    # Now index the TOP of YEAR 01-01 and TOP Of MONTH
    f = "%m-%d-%Y"
    list_bop_str = [None]*datlen
    
    for i in range(0,datlen):
        list_bop_str[i] = list_bop[i].strftime(f)
    
    xtick_NY = []
    xtick_NM = []
    xtick_NY_label = []
    xtick_NM_label = []
    
    for j in range(0,datlen):
        if list_bop_str[j][0:5]=='01-01':
            xtick_NY.append(j)
            xtick_NY_label.append(list_bop_str[j][6:10]) # YYYY
        if list_bop_str[j][3:5]=='01':
            xtick_NM.append(j)
            xtick_NM_label.append(list_bop_str[j][0:2]) # MM
    return xtick_NY, xtick_NY_label, xtick_NM, xtick_NM_label


#---------------------------------------------------------------------------------
# Create Vault for time series of observed time series ... 
#--------------------------------------------------------------------------------




    
def make_hydrologic_plots(VaultData,mlrun,sprun,ml_info):
    # -----------------------------------
    N =  (ml_info.edate - ml_info.sdate).days
    key_intvl = ml_info.keyintvl
    
    sim1 = ml_info.sdate
    sim2 = ml_info.edate
    # ---------------------
    # Grab Data from Vault
    # --------------------
    
    sp_nbs = VaultData.withdraw(kind='nbs', units='cms', intvl=key_intvl, loc='sup', first=sim1, last=sim2).dataVals
    sp_div = VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='oll', first=sim1, last=sim2).dataVals
    sp_qou = VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='smr', first=sim1, last=sim2).dataVals    
    sp_dst = sprun[0]
    sp_dstN = [-x for x in sp_dst]
    
    mh_nbs = VaultData.withdraw(kind='nbs', units='cms', intvl=key_intvl, loc='mhu', first=sim1, last=sim2).dataVals
    mh_qin = VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='smr', first=sim1, last=sim2).dataVals
    mh_div = VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='chi', first=sim1, last=sim2).dataVals
    mh_qou = VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='scr', first=sim1, last=sim2).dataVals
    mh_dst = mlrun[0] # storage piece
    mh_dstN = [-x for x in mh_dst]
    
    N = len(mh_nbs)
    ind = np.arange(N)
        
    sc_nbs = VaultData.withdraw(kind='nbs', units='cms', intvl=key_intvl, loc='stc', first=sim1, last=sim2).dataVals
    sc_qin = VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='scr', first=sim1, last=sim2).dataVals
    sc_div = [0]*N
    sc_qou = VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='det', first=sim1, last=sim2).dataVals
    sc_dst = mlrun[1]
    sc_dstN = [-x for x in sc_dst]

    er_nbs = VaultData.withdraw(kind='nbs', units='cms', intvl=key_intvl, loc='eri', first=sim1, last=sim2).dataVals
    er_qin = VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='det', first=sim1, last=sim2).dataVals
    er_div = VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='wel', first=sim1, last=sim2).dataVals
    er_qou = VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='nia', first=sim1, last=sim2).dataVals
    er_dst = mlrun[2]
    er_dstN = [-x for x in er_dst]
    
    # ----------------------
    # Make Necessary Series for Plotting
    #-----------------------
    sp_LHS = [None]*N
    sp_RHS = [None]*N
    sp_lmr = [None]*N
    sp_error = [None]*N

    for x in range(0,N):
        sp_LHS[x] = sp_nbs[x] + sp_div[x]
    for x in range(0,N):
        sp_RHS[x] = sp_qou[x]
    for x in range(0,N):
        sp_lmr[x] = sp_LHS[x] - sp_RHS[x]
    for x in range(0,N):
        sp_error[x] = sp_lmr[x] - sp_dst[x]    
    
    mh_LHS = [None]*N
    mh_RHS = [None]*N
    mh_lmr = [None]*N
    mh_error = [None]*N

    for x in range(0,N):
        mh_LHS[x] = mh_nbs[x] + mh_qin[x]
    for x in range(0,N):
        mh_RHS[x] = mh_div[x] + mh_qou[x]
    for x in range(0,N):
        mh_lmr[x] = mh_LHS[x] - mh_RHS[x]
    for x in range(0,N):
        mh_error[x] = mh_lmr[x] - mh_dst[x]
    
    sc_LHS = [None]*N
    sc_RHS = [None]*N
    sc_lmr = [None]*N
    sc_error = [None]*N
    
    for x in  range(0,N):
        sc_LHS[x] = sc_nbs[x] + sc_qin[x]
    for x in  range(0,N):
        sc_RHS[x] = sc_div[x] + sc_qou[x]
    for x in  range(0,N):
        sc_lmr[x] = sc_LHS[x] - sc_RHS[x]
    for x in  range(0,N):
        sc_error[x] = sc_lmr[x] - sc_dst[x]
    
    er_LHS = [None]*N
    er_RHS = [None]*N
    er_lmr = [None]*N
    er_error = [None]*N
    
    for x in  range(0,N):
        er_LHS[x] = er_nbs[x] + er_qin[x]
    for x in  range(0,N):
        er_RHS[x] = er_div[x] + er_qou[x]
    for x in  range(0,N):
        er_lmr[x] = er_LHS[x] - er_RHS[x]
    for x in  range(0,N):
        er_error[x] = er_lmr[x] - er_dst[x] 
    #---------------------------------
    # DEFINE SOME LIMITS FOR PLOTTING (rounding to 500 cms)
    #---------------------------------

    spinmaxr = max(mh_LHS)
    spinmax = round(spinmaxr/500)*500
    yscalehu = spinmax/2
    
    mhinmaxr = max(mh_LHS)
    mhinmax = round(mhinmaxr/500)*500
    yscalehm = mhinmax/2
    
    
    scinmaxr = max(sc_LHS)
    scinmax = round(scinmaxr/500)*500
    yscalehs = scinmax/2

    
    erinmaxr = max(er_LHS)
    erinmax = round(erinmaxr/500)*500
    yscalehe = erinmax/2

    #-------------------------------
    # GET THE X TICK INFO
    #--------------------------------
    tickYear = daterangeBOP(VaultData,ml_info)[0]
    tickYear_lab = daterangeBOP(VaultData,ml_info)[1]
    
    tickMonth = daterangeBOP(VaultData,ml_info)[2]
    tickMonth_lab = daterangeBOP(VaultData,ml_info)[3]
    
    #-------------------------------
    # TITLE INFO FOR PLOTS & SAVE INFO
    #--------------------------------
    title_sup = 'Hydrologic Equations: Lake Superior \n NBS + D = Qout + dS/dt [cms]'
    save_sup = "data-out-plots/" +  key_intvl + "/Hydrologic_SUP.jpeg"
    
    title_mhu = 'Hydrologic Equations: Lake Michigan-Huron \n NBS + Qin = D + Qout + dS/dt [cms]'
    save_mhu = "data-out-plots/" +  key_intvl + "/Hydrologic_MHU.jpeg"
    
    title_stc = 'Hydrologic Equations: Lake St. Clair \n NBS + Qin =     Qout + dS/dt [cms]'
    save_stc = "data-out-plots/" +  key_intvl + "/Hydrologic_STC.jpeg"
    
    title_eri = 'Hydrologic Equations: Lake Erie \n NBS + Qin = D + Qout + dS/dt [cms]'
    save_eri = "data-out-plots/" +  key_intvl + "/Hydrologic_ERI.jpeg"
    #-------------------------
    # PLOTS: CONTINUITY CHECK (4)
    #-------------------------
    fig, axes = plt.subplots(4, 1, figsize=(8.5,11), sharex=True,gridspec_kw={'hspace':0})
    ax0, ax1, ax2, ax3 = axes.flatten()
    ax0.bar(ind, sp_nbs, width = 1.0, alpha=0.4, color='grey', label='NBS')
    ax0.bar(ind, sp_div, width = 1.0, alpha=0.4, color='b', label='Div', bottom = sp_nbs)
    ax0.step(sp_LHS,'k-',where='mid',linewidth=0.5,label='Net Inflow')
    ax1.bar(ind, sp_qou, width = 1.0, alpha=0.4, color='m', label='Qou')
    #ax1.bar(ind, mh_qou, width = 1.0, alpha=0.4, color='m', label='Qou', bottom = mh_div)
    ax1.step(sp_RHS,'k-',where='mid',linewidth=0.5,label='Net Outflow')
    ax2.step(sp_lmr,'k-',where='mid',linewidth=0.5,label='Inflow - Outflow')
    ax2.bar(ind, sp_dstN, width = 1.0, alpha=0.4, color='c',label='Change in Storage',bottom=sp_lmr)
    ax2.plot([0.0000]*N,'k-',linewidth=0.5,linestyle='dotted')
    ax3.bar(ind, sp_error, width = 1.0, alpha=0.4, color='red',label='Error')
    ax3.plot([0.0000]*N,'k-',linewidth=0.5,linestyle='dotted')
    ax0.legend(); ax1.legend(); ax2.legend(); ax3.legend()
    ax0.set_ylabel('Inflow')
    ax1.set_ylabel('Outflow')
    ax2.set_ylabel('I - O - dS/dt')
    ax3.set_ylabel('Continuity Error')
    ax0.set_ylim([0,spinmax])
    ax1.set_ylim([0,spinmax])
    ax1.invert_yaxis()
    ax3.set_ylim([-yscalehu,yscalehu])
    ax0.yaxis.tick_right();ax1.yaxis.tick_right()
    ax2.yaxis.tick_right();ax3.yaxis.tick_right()
    ax0.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax3.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax0.tick_params(labelsize=8);ax1.tick_params(labelsize=8);
    ax2.tick_params(labelsize=8);ax3.tick_params(labelsize=8)
    ax3.xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax3.xaxis.set_minor_locator(ticker.FixedLocator(tickMonth))
    ax3.xaxis.set_ticklabels(tickYear_lab)
    ax0.set_title(title_sup)
    fig.savefig(save_sup,dpi=(200))
    
    fig, axes = plt.subplots(4, 1, figsize=(8.5,11), sharex=True,gridspec_kw={'hspace':0})
    ax0, ax1, ax2, ax3 = axes.flatten()
    ax0.bar(ind, mh_qin,  width = 1.0, alpha=0.4, color='b', label='Qin')
    ax0.bar(ind, mh_nbs,  width = 1.0, alpha=0.4, color='grey', label='NBS', bottom = mh_qin)
    ax0.step(mh_LHS,'k-',where='mid',linewidth=0.5,label='Net Inflow')
    ax1.bar(ind, mh_div, width = 1.0, alpha=0.4, color='darkgrey', label='Div')
    ax1.bar(ind, mh_qou, width = 1.0, alpha=0.4, color='m', label='Qou', bottom = mh_div)
    ax1.step(mh_RHS,'k-',where='mid',linewidth=0.5,label='Net Outflow')
    ax2.step(mh_lmr,'k-',where='mid',linewidth=0.5,label='Inflow - Outflow')
    ax2.bar(ind, mh_dstN, width = 1.0, alpha=0.4, color='c',label='Change in Storage',bottom=mh_lmr)
    ax2.plot([0.0000]*N,'k-',linewidth=0.5,linestyle='dotted')
    ax3.bar(ind, mh_error, width = 1.0, alpha=0.4, color='red',label='Error')
    ax3.plot([0.0000]*N,'k-',linewidth=0.5,linestyle='dotted')
    ax0.legend(); ax1.legend(); ax2.legend(); ax3.legend()
    ax0.set_ylabel('Inflow')
    ax1.set_ylabel('Outflow')
    ax2.set_ylabel('I - O - dS/dt')
    ax3.set_ylabel('Continuity Error')
    ax0.set_ylim([0,mhinmax])
    ax1.set_ylim([0,mhinmax])
    ax1.invert_yaxis()
    ax3.set_ylim([-yscalehm,yscalehm])
    ax0.yaxis.tick_right();ax1.yaxis.tick_right()
    ax2.yaxis.tick_right();ax3.yaxis.tick_right()
    ax0.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax3.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax0.tick_params(labelsize=8);ax1.tick_params(labelsize=8);
    ax2.tick_params(labelsize=8);ax3.tick_params(labelsize=8)
    ax3.xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax3.xaxis.set_minor_locator(ticker.FixedLocator(tickMonth))
    ax3.xaxis.set_ticklabels(tickYear_lab)
    ax0.set_title(title_mhu)
    fig.savefig(save_mhu,dpi=(200))

    fig, axes = plt.subplots(4, 1, figsize=(8.5,11), sharex=True,gridspec_kw={'hspace':0})
    ax0, ax1, ax2, ax3 = axes.flatten()
    ax0.bar(ind, sc_qin,  width = 1.0, alpha=0.4, color='b', label='Qin')
    ax0.bar(ind, sc_nbs,  width = 1.0, alpha=0.4, color='grey', label='NBS', bottom = sc_qin)
    ax0.step(sc_LHS,'k-',where='mid',linewidth=0.5,label='Net Inflow')
    ax1.bar(ind, sc_div, width = 1.0, alpha=0.4, color='darkgrey', label='Div')
    ax1.bar(ind, sc_qou, width = 1.0, alpha=0.4, color='m', label='Qou', bottom = sc_div)
    ax1.step(sc_RHS,'k-',where='mid',linewidth=0.5,label='Net Outflow')
    ax2.step(sc_lmr,'k-',where='mid',linewidth=0.5,label='Inflow - Outflow')
    ax2.bar(ind, sc_dstN, width = 1.0, alpha=0.4, color='c',label='Change in Storage',bottom=sc_lmr)
    ax2.plot([0.0000]*N,'k-',linewidth=0.5,linestyle='dotted')
    ax3.bar(ind, sc_error, width = 1.0, alpha=0.4, color='red',label='Error')
    ax3.plot([0.0000]*N,'k-',linewidth=0.5,linestyle='dotted')
    ax0.legend(); ax1.legend(); ax2.legend(); ax3.legend()
    ax0.set_ylabel('Inflow')
    ax1.set_ylabel('Outflow')
    ax2.set_ylabel('I - O - dS/dt')
    ax3.set_ylabel('Continuity Error')
    ax0.set_ylim([0,scinmax])
    ax1.set_ylim([0,scinmax])
    ax1.invert_yaxis()
    ax3.set_ylim([-yscalehs,yscalehs])
    ax0.yaxis.tick_right();ax1.yaxis.tick_right()
    ax2.yaxis.tick_right();ax3.yaxis.tick_right()
    ax0.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax3.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax0.tick_params(labelsize=8);ax1.tick_params(labelsize=8);
    ax2.tick_params(labelsize=8);ax3.tick_params(labelsize=8)
    ax3.xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax3.xaxis.set_minor_locator(ticker.FixedLocator(tickMonth))
    ax3.xaxis.set_ticklabels(tickYear_lab)
    ax0.set_title(title_stc)
    fig.savefig(save_stc,dpi=(200))    
    
    fig, axes = plt.subplots(4, 1, figsize=(8.5,11), sharex=True,gridspec_kw={'hspace':0})
    ax0, ax1, ax2, ax3 = axes.flatten()
    ax0.bar(ind, er_qin,  width = 1.0, alpha=0.4, color='b', label='Qin')
    ax0.bar(ind, er_nbs,  width = 1.0, alpha=0.4, color='grey', label='NBS', bottom = er_qin)
    ax0.step(er_LHS,'k-',where='mid',linewidth=0.5,label='Net Inflow')
    ax1.bar(ind, er_div, width = 1.0, alpha=0.4, color='darkgrey', label='Div')
    ax1.bar(ind, er_qou, width = 1.0, alpha=0.4, color='m', label='Qou', bottom = er_div)
    ax1.step(er_RHS,'k-',where='mid',linewidth=0.5,label='Net Outflow')
    ax2.step(er_lmr,'k-',where='mid',linewidth=0.5,label='Inflow - Outflow')
    ax2.bar(ind, er_dstN, width = 1.0, alpha=0.4, color='c',label='Change in Storage',bottom=er_lmr)
    ax2.plot([0.0000]*N,'k-',linewidth=0.5,linestyle='dotted')
    ax3.bar(ind, er_error, width = 1.0, alpha=0.4, color='red',label='Error')
    ax3.plot([0.0000]*N,'k-',linewidth=0.5,linestyle='dotted')
    ax0.legend(); ax1.legend(); ax2.legend(); ax3.legend()
    ax0.set_ylabel('Inflow')
    ax1.set_ylabel('Outflow')
    ax2.set_ylabel('I - O - dS/dt')
    ax3.set_ylabel('Continuity Error')
    ax0.set_ylim([0,erinmax])
    ax1.set_ylim([0,erinmax])
    ax1.invert_yaxis()
    ax3.set_ylim([-yscalehe,yscalehe])
    ax0.yaxis.tick_right();ax1.yaxis.tick_right()
    ax2.yaxis.tick_right();ax3.yaxis.tick_right()
    ax0.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax3.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax0.tick_params(labelsize=8);ax1.tick_params(labelsize=8);
    ax2.tick_params(labelsize=8);ax3.tick_params(labelsize=8)
    ax3.xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax3.xaxis.set_minor_locator(ticker.FixedLocator(tickMonth))
    ax3.xaxis.set_ticklabels(tickYear_lab)
    ax0.set_title(title_eri)
    fig.savefig(save_eri,dpi=(200))
    
    # ------------------------------
    return 'Mass Plots of Hydrologic Equations - see ' + key_intvl + 'folder'

def make_compare_plots(VaultData,VaultObserved,ml_info):
    #-------------------------------
    # GET TIME STEP INFO
    #---------------------------------------
    key_intvl = ml_info.keyintvl
    labelR = ml_info.keyintvl[0].upper()
    plot_header1 = 'Middle Lakes'
    plot_header2 = str(ml_info.titleI)
    plot_header3 = str(ml_info.titleII)
    plot_header4 = 'St Marys Forcing: data'
    
    sim_1 = ml_info.sdate
    sim_2 = ml_info.edate
    
    sim_len = len(VaultObserved.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='sup',first=sim_1,last=sim_2).dataVals)
    
    #-------------------------------
    # FOUR LISTS FULL OF THE DATA
    #--------------------------------    
    
    lev_dat = [VaultObserved.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='sup',first=sim_1,last=sim_2).dataVals,
               VaultObserved.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='mhu',first=sim_1,last=sim_2).dataVals,
               VaultObserved.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='stc',first=sim_1,last=sim_2).dataVals,
               VaultObserved.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='eri',first=sim_1,last=sim_2).dataVals,
               VaultObserved.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='ont',first=sim_1,last=sim_2).dataVals]
    
    flo_dat = [VaultObserved.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='smr',first=sim_1,last=sim_2).dataVals,
               VaultObserved.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='scr',first=sim_1,last=sim_2).dataVals,
               VaultObserved.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='det',first=sim_1,last=sim_2).dataVals,
               VaultObserved.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='nia',first=sim_1,last=sim_2).dataVals,
               VaultObserved.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='stl',first=sim_1,last=sim_2).dataVals,]
    
    # [183.5]*sim_len, [75.0]*sim_len
    lev_res = [VaultData.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='sup',first=sim_1,last=sim_2).dataVals,
               VaultData.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='mhu',first=sim_1,last=sim_2).dataVals,
               VaultData.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='stc',first=sim_1,last=sim_2).dataVals,
               VaultData.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='eri',first=sim_1,last=sim_2).dataVals,
               VaultData.withdraw(kind='mlv', units='m', intvl=key_intvl, loc='ont',first=sim_1,last=sim_2).dataVals]
    
    # [2500]*sim_len , [7500]*sim_len
    flo_res = [VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='smr',first=sim_1,last=sim_2).dataVals,
               VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='scr',first=sim_1,last=sim_2).dataVals,
               VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='det',first=sim_1,last=sim_2).dataVals,
               VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='nia',first=sim_1,last=sim_2).dataVals,
               VaultData.withdraw(kind='flw', units='cms', intvl=key_intvl, loc='stl',first=sim_1,last=sim_2).dataVals]

    #-------------------------------
    # CALCULATE THE ERRORS (mm) & (cms)
    #--------------------------------  
    N = 5
    
    empty_list = [None]*sim_len
    
    bias_lev = [empty_list]*N
    bias_flo = [empty_list]*N
  
    
    for x in range(0,N):
        list_lev_res = lev_res[x]
        list_lev_dat = lev_dat[x]
        diff_lev = [None]*sim_len
    
        for y in range(0,sim_len):
            # convert from m to cm
            mtomm = 100
            diff_lev[y] = (list_lev_res[y] - list_lev_dat[y])*mtomm
            #print(nametemp)
        bias_lev[x] = diff_lev
            
        
    for x in range(0,N):
        list_flo_res = flo_res[x]
        list_flo_dat = flo_dat[x]
        diff_flo = [None]*sim_len
     
        for y in range(0,sim_len):
            diff_flo[y] = list_flo_res[y] - list_flo_dat[y]
            #print(nametemp)
        bias_flo[x] = diff_flo     
            
    #print(bias_lev)
    #print(bias_flo)
    
    bias_level_max = round(max(max(bias_lev)),-1) + 1
    bias_flow_max = round(max(max(bias_flo)),-1) + 1
    
    bias_level_min = round(min(min(bias_lev)),-1) - 10
    bias_flow_min = round(min(min(bias_flo)),-1) - 10
    
    bias_level_lim = max(abs(bias_level_min),abs(bias_level_max))
    bias_flow_lim = max(abs(bias_flow_min),abs(bias_flow_max))
    #-------------------------------
    # GET THE X TICK INFO
    #--------------------------------
    
    
    tickYear = daterangeBOP(VaultData,ml_info)[0]
    tickYear_lab = daterangeBOP(VaultData,ml_info)[1]
    
    tickMonth = daterangeBOP(VaultData,ml_info)[2]
    tickMonth_lab = daterangeBOP(VaultData,ml_info)[3]
    
    
    #-------------------------------
    # TITLE INFO FOR PLOTS & SAVE INFO
    #--------------------------------
    title_levels = 'Great Lakes Regulation and Routing Model (GLRRM): \n Levels (m) \n Superior: dummy, Middle: Iterative, Ontario: dummy'
    save_levels = "data-out-plots/" +  key_intvl + "/GLRRM_levels_&_observed.jpeg"
    
    title_flows = 'Great Lakes Regulation and Routing Model (GLRRM) \n Flows (cms) \n Superior: dummy, Middle: Iterative, Ontario: dummy'
    save_flows = "data-out-plots/" +  key_intvl + "/GLRRM_flows_&_observed.jpeg"
    
    title_levels_bias = 'Great Lakes Regulation and Routing Model (GLRRM) : \n Level bias (cm)'
    save_levels_bias = "data-out-plots/" +  key_intvl + "/Bias_GLRRM_levels.jpeg"
    
    title_flows_bias = 'Great Lakes Regulation and Routing Model (GLRRM) : \n Flows bias (cms)'
    save_flows_bias = "data-out-plots/" +  key_intvl + "/Bias_GLRRM_flows.jpeg"
    #-------------------------
    # PLOTS: SIMULATED vs OBSERVED LEVELS & FLOWS (2)
    #-------------------------

    #print(tickYear)
    #print(tickYear_lab)
    
    fig, ax = plt.subplots(5,sharex=True,figsize=(8.5,11),gridspec_kw={'hspace':0})
    ax[0].plot(lev_dat[0],linestyle='-', color='black',drawstyle='steps-post',label='observed')
    ax[0].plot(lev_res[0],linestyle='--', color='magenta',drawstyle='steps-post',label='glrrm-dummy')
    ax[1].plot(lev_dat[1],linestyle='-', color='black',drawstyle='steps-post',label='observed')
    ax[1].plot(lev_res[1],linestyle='-', color='c',drawstyle='steps-post',label='glrrm')
    ax[2].plot(lev_dat[2],linestyle='-', color='black',drawstyle='steps-post',label='observed')
    ax[2].plot(lev_res[2],linestyle='-', color='c',drawstyle='steps-post',label='glrrm')
    ax[3].plot(lev_dat[3],linestyle='-', color='black',drawstyle='steps-post',label='observed')
    ax[3].plot(lev_res[3],linestyle='-', color='c',drawstyle='steps-post',label='glrrm')
    ax[4].plot(lev_dat[4],linestyle='-', color='black',drawstyle='steps-post',label='observed')
    ax[4].plot(lev_res[4],linestyle='--', color='lawngreen',drawstyle='steps-post',label='glrrm-dummy')
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    ax[3].legend()
    ax[4].legend()
    ax[0].set_ylabel('superior')
    ax[1].set_ylabel('michigan-huron')
    ax[2].set_ylabel('st.clair')
    ax[3].set_ylabel('erie')
    ax[4].set_ylabel('ontario')
    ax[0].yaxis.tick_right();ax[1].yaxis.tick_right();ax[2].yaxis.tick_right();
    ax[3].yaxis.tick_right();ax[4].yaxis.tick_right()
    ax[0].tick_params(labelsize=8);ax[1].tick_params(labelsize=8);ax[2].tick_params(labelsize=8)
    ax[3].tick_params(labelsize=8);ax[4].tick_params(labelsize=8)
    ax[0].set_title(title_levels)
    ax[4].xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax[4].xaxis.set_minor_locator(ticker.FixedLocator(tickMonth))
    ax[4].xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax[4].xaxis.set_ticklabels(tickYear_lab)
    plt.subplots_adjust(top=0.83)
    plt.show()
    fig.savefig(save_levels,dpi=(200))


    fig, ax = plt.subplots(5,sharex=True,figsize=(8.5,11),gridspec_kw={'hspace':0})
    ax[0].plot(flo_dat[0],linestyle='-', color='black',drawstyle='steps-post',label='observed')
    ax[0].plot(flo_res[0],linestyle='--',color='magenta',drawstyle='steps-post',label='glrrm-dummy')
    ax[1].plot(flo_dat[1],linestyle='-',color='black',drawstyle='steps-post',label='observed')
    ax[1].plot(flo_res[1],linestyle='-',color='cyan',drawstyle='steps-post',label='glrrm')
    ax[2].plot(flo_dat[2],linestyle='-',color='black',drawstyle='steps-post',label='observed')
    ax[2].plot(flo_res[2],linestyle='-',color='cyan',drawstyle='steps-post',label='glrrm')
    ax[3].plot(flo_dat[3],linestyle='-',color='black',drawstyle='steps-post',label='observed')
    ax[3].plot(flo_res[3],linestyle='-',color='cyan',drawstyle='steps-post',label='glrrm')
    ax[4].plot(flo_dat[4],linestyle='-',color='black',drawstyle='steps-post',label='observed')
    ax[4].plot(flo_res[4],linestyle='--',color='lawngreen',drawstyle='steps-post',label='glrrm-dummy')
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    ax[3].legend()
    ax[4].legend()
    ax[0].set_ylabel('st.marys river')
    ax[1].set_ylabel('st.clair river')
    ax[2].set_ylabel('detroit river')
    ax[3].set_ylabel('niagara river')
    ax[4].set_ylabel('st.lawrence river')
    ax[0].yaxis.tick_right();ax[1].yaxis.tick_right();ax[2].yaxis.tick_right();
    ax[3].yaxis.tick_right();ax[4].yaxis.tick_right()
    ax[0].tick_params(labelsize=8);ax[1].tick_params(labelsize=8);ax[2].tick_params(labelsize=8)
    ax[3].tick_params(labelsize=8);ax[4].tick_params(labelsize=8)
    ax[0].set_title(title_flows)
    ax[4].xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax[4].xaxis.set_minor_locator(ticker.FixedLocator(tickMonth))
    ax[4].xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax[4].xaxis.set_ticklabels(tickYear_lab)
    plt.subplots_adjust(top=0.83)
    plt.show()
    fig.savefig(save_flows,dpi=(200))
    
    #-------------------------
    # PLOTS: GLRRM BIAS (mm & csm) (2)
    #-------------------------
    ind = np.arange(sim_len)
        
    fig, ax = plt.subplots(5,sharex=True,figsize=(8.5,11),gridspec_kw={'hspace':0})
    ax[0].bar(ind,bias_lev[0],width=1.0,color='blue',alpha=0.4,label='glrrm-bias')
    ax[1].bar(ind,bias_lev[1],width=1.0,color='blue',alpha=0.4,label='glrrm-bias')
    ax[2].bar(ind,bias_lev[2],width=1.0,color='blue',alpha=0.4,label='glrrm-bias')
    ax[3].bar(ind,bias_lev[3],width=1.0,color='blue',alpha=0.4,label='glrrm-bias')
    ax[4].bar(ind,bias_lev[4],width=1.0,color='blue',alpha=0.4,label='glrrm-bias')
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    ax[3].legend()
    ax[4].legend()
    ax[0].set_ylabel('superior')
    ax[1].set_ylabel('michigan-huron')
    ax[2].set_ylabel('st.clair')
    ax[3].set_ylabel('erie')
    ax[4].set_ylabel('ontario')
    ax[0].yaxis.tick_right();ax[1].yaxis.tick_right();ax[2].yaxis.tick_right();
    ax[3].yaxis.tick_right();ax[4].yaxis.tick_right()
    #ax[0].set_ylim([-bias_level_lim,bias_level_lim])
    #ax[1].set_ylim([-bias_level_lim,bias_level_lim])
    #ax[2].set_ylim([-bias_level_lim,bias_level_lim])
    #ax[3].set_ylim([-bias_level_lim,bias_level_lim])
    #ax[4].set_ylim([-bias_level_lim,bias_level_lim])
    ax[0].tick_params(labelsize=8);ax[1].tick_params(labelsize=8);ax[2].tick_params(labelsize=8)
    ax[3].tick_params(labelsize=8);ax[4].tick_params(labelsize=8)
    ax[0].set_title(title_levels_bias)
    ax[4].xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax[4].xaxis.set_minor_locator(ticker.FixedLocator(tickMonth))
    ax[4].xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax[4].xaxis.set_ticklabels(tickYear_lab)
    plt.subplots_adjust(top=0.83)
    plt.show()
    fig.savefig(save_levels_bias,dpi=(200))

    fig, ax = plt.subplots(5,sharex=True,figsize=(8.5,11),gridspec_kw={'hspace':0})
    ax[0].bar(ind,bias_flo[0],width=1.0,color='red',alpha=0.4,label='glrrm-bias')
    ax[1].bar(ind,bias_flo[1],width=1.0,color='red',alpha=0.4,label='glrrm-bias')
    ax[2].bar(ind,bias_flo[2],width=1.0,color='red',alpha=0.4,label='glrrm-bias')
    ax[3].bar(ind,bias_flo[3],width=1.0,color='red',alpha=0.4,label='glrrm-bias')
    ax[4].bar(ind,bias_flo[4],width=1.0,color='red',alpha=0.4,label='glrrm-bias')
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    ax[3].legend()
    ax[4].legend()
    ax[0].set_ylabel('superior')
    ax[1].set_ylabel('michigan-huron')
    ax[2].set_ylabel('st.clair')
    ax[3].set_ylabel('erie')
    ax[4].set_ylabel('ontario')
    ax[0].yaxis.tick_right();ax[1].yaxis.tick_right();ax[2].yaxis.tick_right();
    ax[3].yaxis.tick_right();ax[4].yaxis.tick_right()
    #ax[0].set_ylim([-bias_flow_lim,bias_flow_lim])
    #ax[1].set_ylim([-bias_flow_lim,bias_flow_lim])
    #ax[2].set_ylim([-bias_flow_lim,bias_flow_lim])
    #ax[3].set_ylim([-bias_flow_lim,bias_flow_lim])
    #ax[4].set_ylim([-bias_flow_lim,bias_flow_lim])
    ax[0].tick_params(labelsize=8);ax[1].tick_params(labelsize=8);ax[2].tick_params(labelsize=8)
    ax[3].tick_params(labelsize=8);ax[4].tick_params(labelsize=8)
    ax[0].set_title(title_flows_bias)
    ax[4].xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax[4].xaxis.set_minor_locator(ticker.FixedLocator(tickMonth))
    ax[4].xaxis.set_major_locator(ticker.FixedLocator(tickYear))
    ax[4].xaxis.set_ticklabels(tickYear_lab)
    plt.subplots_adjust(top=0.83)
    plt.show()
    fig.savefig(save_flows_bias,dpi=(200))


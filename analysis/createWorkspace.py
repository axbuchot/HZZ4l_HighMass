import ROOT
import os,sys,subprocess
from math import trunc

from paths import path

# Considering or not decimals in bin boundaries
decimal = {
'mass4l': False,
'mass4l_zzfloating': False,
'njets_pt30_eta4p7': False,
'pT4l': False,
'rapidity4l': True,
'costhetaZ1': True,
'costhetaZ2': True,
'phi': True,
'phistar': True,
'costhetastar': True,
'massZ1': False,
'massZ2': False,
'pTj1': False,
'pTHj': False,
'mHj': False,
'pTj2': False,
'mjj': False,
'absdetajj': True,
'dphijj': True,
'pTHjj': False,
'TCjmax': False,
'TBjmax': False,
'D0m': True,
'Dcp': True,
'D0hp': True,
'Dint': True,
'DL1': True,
'DL1Zg': True,
'rapidity4l_pT4l': True,
'njets_pt30_eta4p7 vs pT4l': False,
'pTj1_pTj2': False,
'pT4l_pTHj': False,
'massZ1_massZ2': False,
'TCjmax_pT4l': False
}


sys.path.append('../../inputs/')
sys.path.append('../../templates/')

def createWorkspace(obsName, channel, nBins, obsBin, observableBins, usecfactor, addfakeH, modelName, physicalModel, year, JES, doubleDiff, lowerBound, upperBound, rawObsName):
    print '\n'
    print 'Creating WorkSpace', year

    if not doubleDiff:
        obsBin_low = observableBins[obsBin]
        obsBin_high = observableBins[obsBin+1]
        obs_bin_lowest = observableBins[0]
        obs_bin_highest = observableBins[len(observableBins)-1]
    elif doubleDiff:
        obsBin_low = observableBins[obsBin][0]
        obsBin_high = observableBins[obsBin][1]
        obs_bin_lowest = min(x[0] for x in observableBins.values())
        obs_bin_highest = max(x[1] for x in observableBins.values())
        #second variables
        obsBin_2nd_low = observableBins[obsBin][2]
        obsBin_2nd_high = observableBins[obsBin][3]
        obs_bin_2nd_lowest = min(x[2] for x in observableBins.values())
        obs_bin_2nd_highest = max(x[3] for x in observableBins.values())

    recobin = "recobin"+str(obsBin)
    print recobin
    JES = False
    doJES = False

    # Load some libraries
    ROOT.gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
    ROOT.gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so") #Print 0 in case of succesfull loading
    ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include")
    ROOT.gSystem.AddIncludePath("-Iinclude/")

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    # from inputs_sig imports some coefficients
    if (usecfactor):
        _temp = __import__('inputs_sig_'+obsName+'_'+year, globals(), locals(), ['cfactor','inc_wrongfrac','binfrac_wrongfrac','inc_outfrac','binfrac_outfrac'], -1)
        cfactor = _temp.cfactor
        inc_outfrac = _temp.inc_outfrac
        binfrac_outfrac = _temp.binfrac_wrongfrac
    else:
        _temp = __import__('inputs_sig_'+obsName+'_'+year, globals(), locals(), ['acc','eff','inc_wrongfrac','binfrac_wrongfrac','outinratio','lambdajesup','lambdajesdn'], -1)
        acc = _temp.acc
        eff = _temp.eff
        outinratio = _temp.outinratio
    lambdajesup = _temp.lambdajesup
    lambdajesdn = _temp.lambdajesdn
    inc_wrongfrac = _temp.inc_wrongfrac
    binfrac_wrongfrac = _temp.binfrac_wrongfrac
    #number_fake = _temp.number_fake

    # import h4l highmass br
    _temp = __import__('higgs_highmass_br_13TeV', globals(), locals(), ['higgs_highmass','higgs4l_br'], -1)
    higgs_highmass = _temp.higgs_highmass
    higgs4l_br = _temp.higgs4l_br

    _temp = __import__('inputs_bkg_'+obsName+'_'+year, globals(), locals(), ['fractionsBackground'], -1)
    fractionsBackground = _temp.fractionsBackground

    _temp = __import__('observables', globals(), locals(), ['observables'], -1)
    observables = _temp.observables

    if (not obsName=="mass4l"):
        obsName_help = observables[rawObsName]['obs_reco']
        if doubleDiff: obsName_2nd_help = observables[rawObsName]['obs_reco_2nd']

        observable = ROOT.RooRealVar(obsName_help,obsName_help,float(obs_bin_lowest),float(obs_bin_highest))
        if doubleDiff: observable_2nd = ROOT.RooRealVar(obsName_2nd_help,obsName_2nd_help,float(obs_bin_2nd_lowest),float(obs_bin_2nd_highest)) #ATdouble

        observable.Print()
        if doubleDiff: observable_2nd.Print()

    # Parameters of doubleCB signal
    m = ROOT.RooRealVar("CMS_zz4l_mass", "CMS_zz4l_mass", lowerBound, upperBound)
    #MH = ROOT.RooRealVar("MH","MH",125.7,109.55,1000.05)
    MH = ROOT.RooRealVar("MH","MH", 125.38, lowerBound, upperBound)
    if(year == '2018'):
        CMS_zz4l_mean_m_sig_2018 = ROOT.RooRealVar("CMS_zz4l_mean_m_sig_2018","CMS_zz4l_mean_m_sig_2018",-10,10)
        CMS_zz4l_mean_e_sig_2018 = ROOT.RooRealVar("CMS_zz4l_mean_e_sig_2018","CMS_zz4l_mean_e_sig_2018",-10,10)
        CMS_zz4l_sigma_m_sig_2018 = ROOT.RooRealVar("CMS_zz4l_sigma_m_sig_2018","CMS_zz4l_sigma_m_sig_2018",-10,10)
        CMS_zz4l_sigma_e_sig_2018 = ROOT.RooRealVar("CMS_zz4l_sigma_e_sig_2018","CMS_zz4l_sigma_e_sig_2018",-10,10)
        lumi = ROOT.RooRealVar("lumi_132018","lumi_132018", 59.7)
        if (channel=='2e2mu'):
            CMS_zz4l_mean_m_err_3_2018 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_3_2018","CMS_zz4l_mean_m_err_3_2018",0.0004,0.0004,0.0004)
            CMS_zz4l_mean_e_err_3_2018 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_3_2018","CMS_zz4l_mean_e_err_3_2018",0.003,0.003,0.003)
            CMS_zz4l_n_sig_3_2018 = ROOT.RooRealVar("CMS_zz4l_n_sig_3_2018","CMS_zz4l_n_sig_3_2018",-10,10)
            CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2018", "(124.260469656+(0.995095874123)*(@0-125)) + (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2018,CMS_zz4l_mean_e_sig_2018,CMS_zz4l_mean_m_err_3_2018,CMS_zz4l_mean_e_err_3_2018))
            CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2018","CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2018","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2018))
            CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2018", "(1.55330758963+(0.00797274642218)*(@0-125))*(TMath::Sqrt((1+@1)*(1+@2)))",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2018,CMS_zz4l_sigma_e_sig_2018))
            CMS_zz4l_alpha_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2018","(0.947414158515+(0)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2018","(3.33147279858+(-0.0438375854704)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_3_2018))
            CMS_zz4l_alpha2_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2018","(1.52497361611+(0)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_3_centralValue_2e2murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2018","CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2018","(5.20522265056+(0)*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB("trueH", "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2018, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2018, CMS_zz4l_alpha_3_centralValue_2e2murecobin2018, CMS_zz4l_n_3_centralValue_2e2murecobin2018, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2018, CMS_zz4l_n2_3_centralValue_2e2murecobin2018)
        if (channel=='4e'):
            CMS_zz4l_mean_e_err_2_2018 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_2_2018","CMS_zz4l_mean_e_err_2_2018",0.003,0.003,0.003)
            CMS_zz4l_n_sig_2_2018 = ROOT.RooRealVar("CMS_zz4l_n_sig_2_2018","CMS_zz4l_n_sig_2_2018",-10,10)
            CMS_zz4l_mean_sig_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2018","CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2018", "(123.5844824+(0.985478630993)*(@0-125)) + @0*@1*@2", ROOT.RooArgList(MH,CMS_zz4l_mean_e_sig_2018,CMS_zz4l_mean_e_err_2_2018))
            CMS_zz4l_mean_sig_NoConv_2_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_2_4erecobin2018","CMS_zz4l_mean_sig_NoConv_2_4erecobin2018","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_2_centralValue_4erecobin2018))
            CMS_zz4l_sigma_sig_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2018","CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2018", "(2.06515102908+(0.0170917403402)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_e_sig_2018))
            CMS_zz4l_alpha_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2018","CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2018","(0.948100247167+(0)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n_2_centralValue_4e" + recobin + "2018","CMS_zz4l_n_2_centralValue_4e" + recobin + "2018","(4.50639853892+(0)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_2_2018))
            CMS_zz4l_alpha2_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2018","CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2018","(1.50095152675+(0)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_2_centralValue_4erecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n2_2_centralValue_4e" + recobin + "2018","CMS_zz4l_n2_2_centralValue_4e" + recobin + "2018","(8.41693578742+(0.219719825966)*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB("trueH", "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2018, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2018, CMS_zz4l_alpha_2_centralValue_4erecobin2018, CMS_zz4l_n_2_centralValue_4erecobin2018, CMS_zz4l_alpha2_2_centralValue_4erecobin2018, CMS_zz4l_n2_2_centralValue_4erecobin2018)
        if (channel=='4mu'):
            CMS_zz4l_mean_m_err_1_2018 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_1_2018","CMS_zz4l_mean_m_err_1_2018",0.0004,0.0004,0.0004)
            CMS_zz4l_n_sig_1_2018 = ROOT.RooRealVar("CMS_zz4l_n_sig_1_2018","CMS_zz4l_n_sig_1_2018",-10,10)
            CMS_zz4l_mean_sig_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2018", "(124.820536957+(0.999619883119)*(@0-125)) + @0*@1*@2",ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2018,CMS_zz4l_mean_m_err_1_2018))
            CMS_zz4l_mean_sig_NoConv_1_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_1_4murecobin2018","CMS_zz4l_mean_sig_NoConv_1_4murecobin2018","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_1_centralValue_4murecobin2018))
            CMS_zz4l_sigma_sig_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2018", "(1.09001384743+(0.00899911411679)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2018))
            CMS_zz4l_alpha_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2018","(1.23329827124+(0)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_n_1_centralValue_4mu" + recobin + "2018","(2.04575884495+(0)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_1_2018))
            CMS_zz4l_alpha2_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2018","(1.84386824883+(0)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_1_centralValue_4murecobin2018 = ROOT.RooFormulaVar("CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2018","CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2018","(2.98483993137+(0)*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB("trueH", "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2018, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2018, CMS_zz4l_alpha_1_centralValue_4murecobin2018, CMS_zz4l_n_1_centralValue_4murecobin2018, CMS_zz4l_alpha2_1_centralValue_4murecobin2018, CMS_zz4l_n2_1_centralValue_4murecobin2018)

        # Wrong signal combination events
        if (channel=='4mu'):
            p1_12018 = ROOT.RooRealVar("CMS_fakeH_p1_12018","p1_12018",165.0, 145.0, 185.0)
            p3_12018 = ROOT.RooRealVar("CMS_fakeH_p3_12018","p3_12018",89.0, 84.0,94.0)
            p2_12018 = ROOT.RooFormulaVar("CMS_fakeH_p2_12018","p2_12018","0.72*@0-@1",ROOT.RooArgList(p1_12018,p3_12018))
            fakeH = ROOT.RooLandau("fakeH", "landau", m, p1_12018, p2_12018)
        if (channel=='4e'):
            p1_22018 = ROOT.RooRealVar("CMS_fakeH_p1_22018","p1_22018",165.0, 145.0, 185.0)
            p3_22018 = ROOT.RooRealVar("CMS_fakeH_p3_22018","p3_22018",89.0, 84.0,94.0)
            p2_22018 = ROOT.RooFormulaVar("CMS_fakeH_p2_22018","p2_22018","0.72*@0-@1",ROOT.RooArgList(p1_22018,p3_22018))
            fakeH = ROOT.RooLandau("fakeH", "landau", m, p1_22018, p2_22018)
        if (channel=='2e2mu'):
            p1_32018 = ROOT.RooRealVar("CMS_fakeH_p1_32018","p1_32018",165.0, 145.0, 185.0)
            p3_32018 = ROOT.RooRealVar("CMS_fakeH_p3_32018","p3_32018",89.0, 84.0,94.0)
            p2_32018 = ROOT.RooFormulaVar("CMS_fakeH_p2_32018","p2_32018","0.72*@0-@1",ROOT.RooArgList(p1_32018,p3_32018))
            fakeH = ROOT.RooLandau("fakeH", "landau", m, p1_32018, p2_32018)

    elif(year == '2017'):
        CMS_zz4l_mean_m_sig_2017 = ROOT.RooRealVar("CMS_zz4l_mean_m_sig_2017","CMS_zz4l_mean_m_sig_2017",-10,10)
        CMS_zz4l_mean_e_sig_2017 = ROOT.RooRealVar("CMS_zz4l_mean_e_sig_2017","CMS_zz4l_mean_e_sig_2017",-10,10)
        CMS_zz4l_sigma_m_sig_2017 = ROOT.RooRealVar("CMS_zz4l_sigma_m_sig_2017","CMS_zz4l_sigma_m_sig_2017",-10,10)
        CMS_zz4l_sigma_e_sig_2017 = ROOT.RooRealVar("CMS_zz4l_sigma_e_sig_2017","CMS_zz4l_sigma_e_sig_2017",-10,10)
        lumi = ROOT.RooRealVar("lumi_132017","lumi_132017", 41.5)
        if (channel=='2e2mu'):
            CMS_zz4l_mean_m_err_3_2017 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_3_2017","CMS_zz4l_mean_m_err_3_2017",0.0004,0.0004,0.0004)
            CMS_zz4l_mean_e_err_3_2017 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_3_2017","CMS_zz4l_mean_e_err_3_2017",0.003,0.003,0.003)
            CMS_zz4l_n_sig_3_2017 = ROOT.RooRealVar("CMS_zz4l_n_sig_3_2017","CMS_zz4l_n_sig_3_2017",-10,10)
            CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2017", "(124.524+(0.00248708+1)*(@0-125)) + (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2017,CMS_zz4l_mean_e_sig_2017,CMS_zz4l_mean_m_err_3_2017,CMS_zz4l_mean_e_err_3_2017))
            CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2017","CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2017","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2017))
            CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2017", "(1.77228+(0.00526163)*(@0-125))*(TMath::Sqrt((1+@1)*(1+@2)))",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2017,CMS_zz4l_sigma_e_sig_2017))
            CMS_zz4l_alpha_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2017","(0.967963+(-0.0047248)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2017","(3.69774+(0)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_3_2017))
            CMS_zz4l_alpha2_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2017","(1.51606+(-0.000272186)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_3_centralValue_2e2murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2017","CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2017","(6.01048+(0)*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB("trueH", "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2017, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2017, CMS_zz4l_alpha_3_centralValue_2e2murecobin2017, CMS_zz4l_n_3_centralValue_2e2murecobin2017, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2017, CMS_zz4l_n2_3_centralValue_2e2murecobin2017)
        if (channel=='4e'):
            CMS_zz4l_mean_e_err_2_2017 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_2_2017","CMS_zz4l_mean_e_err_2_2017",0.003,0.003,0.003)
            CMS_zz4l_n_sig_2_2017 = ROOT.RooRealVar("CMS_zz4l_n_sig_2_2017","CMS_zz4l_n_sig_2_2017",-10,10)
            CMS_zz4l_mean_sig_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2017","CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2017", "(124.1+(-0.00262293+1)*(@0-125)) + @0*@1*@2", ROOT.RooArgList(MH,CMS_zz4l_mean_e_sig_2017,CMS_zz4l_mean_e_err_2_2017))
            CMS_zz4l_mean_sig_NoConv_2_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_2_4erecobin2017","CMS_zz4l_mean_sig_NoConv_2_4erecobin2017","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_2_centralValue_4erecobin2017))
            CMS_zz4l_sigma_sig_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2017","CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2017", "(2.38283+(0.0155)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_e_sig_2017))
            CMS_zz4l_alpha_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2017","CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2017","(0.972669+(-0.00597402)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n_2_centralValue_4e" + recobin + "2017","CMS_zz4l_n_2_centralValue_4e" + recobin + "2017","(5.05142+(0)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_2_2017))
            CMS_zz4l_alpha2_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2017","CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2017","(1.62625+(0.0121146)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_2_centralValue_4erecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n2_2_centralValue_4e" + recobin + "2017","CMS_zz4l_n2_2_centralValue_4e" + recobin + "2017","(6.30057+(0)*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB("trueH", "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2017, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2017, CMS_zz4l_alpha_2_centralValue_4erecobin2017, CMS_zz4l_n_2_centralValue_4erecobin2017, CMS_zz4l_alpha2_2_centralValue_4erecobin2017, CMS_zz4l_n2_2_centralValue_4erecobin2017)
        if (channel=='4mu'):
            CMS_zz4l_mean_m_err_1_2017 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_1_2017","CMS_zz4l_mean_m_err_1_2017",0.0004,0.0004,0.0004)
            CMS_zz4l_n_sig_1_2017 = ROOT.RooRealVar("CMS_zz4l_n_sig_1_2017","CMS_zz4l_n_sig_1_2017",-10,10)
            CMS_zz4l_mean_sig_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2017", "(124.82+(-0.000560694+1)*(@0-125)) + @0*@1*@2",ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2017,CMS_zz4l_mean_m_err_1_2017))
            CMS_zz4l_mean_sig_NoConv_1_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_1_4murecobin2017","CMS_zz4l_mean_sig_NoConv_1_4murecobin2017","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_1_centralValue_4murecobin2017))
            CMS_zz4l_sigma_sig_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2017", "(1.16647+(0.0124833)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2017))
            CMS_zz4l_alpha_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2017","(1.22997+(0.00256332)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_n_1_centralValue_4mu" + recobin + "2017","(2.07185+(0)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_1_2017))
            CMS_zz4l_alpha2_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2017","(1.92338+(0.0109082)*(@0-125))",ROOT.RooArgList(MH))
            CMS_zz4l_n2_1_centralValue_4murecobin2017 = ROOT.RooFormulaVar("CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2017","CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2017","(2.90336+(0)*(@0-125))",ROOT.RooArgList(MH))
            trueH = ROOT.RooDoubleCB("trueH", "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2017, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2017, CMS_zz4l_alpha_1_centralValue_4murecobin2017, CMS_zz4l_n_1_centralValue_4murecobin2017, CMS_zz4l_alpha2_1_centralValue_4murecobin2017, CMS_zz4l_n2_1_centralValue_4murecobin2017)

        # Wrong signal combination events
        if (channel=='4mu'):
            p1_12017 = ROOT.RooRealVar("CMS_fakeH_p1_12017","p1_12017",165.0, 145.0, 185.0)
            p3_12017 = ROOT.RooRealVar("CMS_fakeH_p3_12017","p3_12017",89.0, 84.0,94.0)
            p2_12017 = ROOT.RooFormulaVar("CMS_fakeH_p2_12017","p2_12017","0.72*@0-@1",ROOT.RooArgList(p1_12017,p3_12017))
            fakeH = ROOT.RooLandau("fakeH", "landau", m, p1_12017, p2_12017)
        if (channel=='4e'):
            p1_22017 = ROOT.RooRealVar("CMS_fakeH_p1_22017","p1_22017",165.0, 145.0, 185.0)
            p3_22017 = ROOT.RooRealVar("CMS_fakeH_p3_22017","p3_22017",89.0, 84.0,94.0)
            p2_22017 = ROOT.RooFormulaVar("CMS_fakeH_p2_22017","p2_22017","0.72*@0-@1",ROOT.RooArgList(p1_22017,p3_22017))
            fakeH = ROOT.RooLandau("fakeH", "landau", m, p1_22017, p2_22017)
        if (channel=='2e2mu'):
            p1_32017 = ROOT.RooRealVar("CMS_fakeH_p1_32017","p1_32017",165.0, 145.0, 185.0)
            p3_32017 = ROOT.RooRealVar("CMS_fakeH_p3_32017","p3_32017",89.0, 84.0,94.0)
            p2_32017 = ROOT.RooFormulaVar("CMS_fakeH_p2_32017","p2_32017","0.72*@0-@1",ROOT.RooArgList(p1_32017,p3_32017))
            fakeH = ROOT.RooLandau("fakeH", "landau", m, p1_32017, p2_32017)

    elif(year == '2016'):
          CMS_zz4l_mean_m_sig_2016 = ROOT.RooRealVar("CMS_zz4l_mean_m_sig_2016","CMS_zz4l_mean_m_sig_2016",-10,10)
          CMS_zz4l_mean_e_sig_2016 = ROOT.RooRealVar("CMS_zz4l_mean_e_sig_2016","CMS_zz4l_mean_e_sig_2016",-10,10)
          CMS_zz4l_sigma_m_sig_2016 = ROOT.RooRealVar("CMS_zz4l_sigma_m_sig_2016","CMS_zz4l_sigma_m_sig_2016",-10,10)
          CMS_zz4l_sigma_e_sig_2016 = ROOT.RooRealVar("CMS_zz4l_sigma_e_sig_2016","CMS_zz4l_sigma_e_sig_2016",-10,10)
          lumi = ROOT.RooRealVar("lumi_132016","lumi_132016", 35.9)
          if (channel=='2e2mu'):
              CMS_zz4l_mean_m_err_3_2016 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_3_2016","CMS_zz4l_mean_m_err_3_2016",0.0004,0.0004,0.0004)
              CMS_zz4l_mean_e_err_3_2016 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_3_2016","CMS_zz4l_mean_e_err_3_2016",0.003,0.003,0.003)
              CMS_zz4l_n_sig_3_2016 = ROOT.RooRealVar("CMS_zz4l_n_sig_3_2016","CMS_zz4l_n_sig_3_2016",-10,10)
              CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_mean_sig_3_centralValue_2e2mu" + recobin + "2016", "(124.539+(-0.00679774+1)*(@0-125)) + (@0*@1*@3 + @0*@2*@4)/2", ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2016,CMS_zz4l_mean_e_sig_2016,CMS_zz4l_mean_m_err_3_2016,CMS_zz4l_mean_e_err_3_2016))
              CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2016","CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2016","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2016))
              CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_sigma_sig_3_centralValue_2e2mu" + recobin + "2016", "(1.64632+(0.016435)*(@0-125))*(TMath::Sqrt((1+@1)*(1+@2)))",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2016,CMS_zz4l_sigma_e_sig_2016))
              CMS_zz4l_alpha_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_alpha_3_centralValue_2e2mu" + recobin + "2016","(0.905389+(0.0029819)*(@0-125))",ROOT.RooArgList(MH))
              CMS_zz4l_n_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_n_3_centralValue_2e2mu" + recobin + "2016","(3.90164+(0)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_3_2016))
              CMS_zz4l_alpha2_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_alpha2_3_centralValue_2e2mu" + recobin + "2016","(1.5737+(0.00776476)*(@0-125))",ROOT.RooArgList(MH))
              CMS_zz4l_n2_3_centralValue_2e2murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2016","CMS_zz4l_n2_3_centralValue_2e2mu" + recobin + "2016","(4.33416+(0)*(@0-125))",ROOT.RooArgList(MH))
              trueH = ROOT.RooDoubleCB("trueH", "DoubleCB", m, CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2016, CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2016, CMS_zz4l_alpha_3_centralValue_2e2murecobin2016, CMS_zz4l_n_3_centralValue_2e2murecobin2016, CMS_zz4l_alpha2_3_centralValue_2e2murecobin2016, CMS_zz4l_n2_3_centralValue_2e2murecobin2016)
          if (channel=='4e'):
              CMS_zz4l_mean_e_err_2_2016 = ROOT.RooRealVar("CMS_zz4l_mean_e_err_2_2016","CMS_zz4l_mean_e_err_2_2016",0.003,0.003,0.003)
              CMS_zz4l_n_sig_2_2016 = ROOT.RooRealVar("CMS_zz4l_n_sig_2_2016","CMS_zz4l_n_sig_2_2016",-10,10)
              CMS_zz4l_mean_sig_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2016","CMS_zz4l_mean_sig_2_centralValue_4e" + recobin + "2016", "(124.194+(-0.0123934+1)*(@0-125)) + @0*@1*@2", ROOT.RooArgList(MH,CMS_zz4l_mean_e_sig_2016,CMS_zz4l_mean_e_err_2_2016))
              CMS_zz4l_mean_sig_NoConv_2_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_2_4erecobin2016","CMS_zz4l_mean_sig_NoConv_2_4erecobin2016","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_2_centralValue_4erecobin2016))
              CMS_zz4l_sigma_sig_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2016","CMS_zz4l_sigma_sig_2_centralValue_4e" + recobin + "2016", "(2.09076+(0.0153247)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_e_sig_2016))
              CMS_zz4l_alpha_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2016","CMS_zz4l_alpha_2_centralValue_4e" + recobin + "2016","(0.778691+(-0.00177387)*(@0-125))",ROOT.RooArgList(MH))
              CMS_zz4l_n_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n_2_centralValue_4e" + recobin + "2016","CMS_zz4l_n_2_centralValue_4e" + recobin + "2016","(6.85936+(0)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_2_2016))
              CMS_zz4l_alpha2_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2016","CMS_zz4l_alpha2_2_centralValue_4e" + recobin + "2016","(1.47389+(0.00503384)*(@0-125))",ROOT.RooArgList(MH))
              CMS_zz4l_n2_2_centralValue_4erecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n2_2_centralValue_4e" + recobin + "2016","CMS_zz4l_n2_2_centralValue_4e" + recobin + "2016","(7.24158+(0)*(@0-125))",ROOT.RooArgList(MH))
              trueH = ROOT.RooDoubleCB("trueH", "DoubleCB", m, CMS_zz4l_mean_sig_2_centralValue_4erecobin2016, CMS_zz4l_sigma_sig_2_centralValue_4erecobin2016, CMS_zz4l_alpha_2_centralValue_4erecobin2016, CMS_zz4l_n_2_centralValue_4erecobin2016, CMS_zz4l_alpha2_2_centralValue_4erecobin2016, CMS_zz4l_n2_2_centralValue_4erecobin2016)
          if (channel=='4mu'):
              CMS_zz4l_mean_m_err_1_2016 = ROOT.RooRealVar("CMS_zz4l_mean_m_err_1_2016","CMS_zz4l_mean_m_err_1_2016",0.0004,0.0004,0.0004)
              CMS_zz4l_n_sig_1_2016 = ROOT.RooRealVar("CMS_zz4l_n_sig_1_2016","CMS_zz4l_n_sig_1_2016",-10,10)
              CMS_zz4l_mean_sig_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_mean_sig_1_centralValue_4mu" + recobin + "2016", "(124.801+(-0.00230642+1)*(@0-125)) + @0*@1*@2",ROOT.RooArgList(MH,CMS_zz4l_mean_m_sig_2016,CMS_zz4l_mean_m_err_1_2016))
              CMS_zz4l_mean_sig_NoConv_1_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_1_4murecobin2016","CMS_zz4l_mean_sig_NoConv_1_4murecobin2016","@0",ROOT.RooArgList(CMS_zz4l_mean_sig_1_centralValue_4murecobin2016))
              CMS_zz4l_sigma_sig_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_sigma_sig_1_centralValue_4mu" + recobin + "2016", "(1.20385+(0.00862539)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_sigma_m_sig_2016))
              CMS_zz4l_alpha_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_alpha_1_centralValue_4mu" + recobin + "2016","(1.29006+(-0.0040219)*(@0-125))",ROOT.RooArgList(MH))
              CMS_zz4l_n_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_n_1_centralValue_4mu" + recobin + "2016","(2.1216+(0)*(@0-125))*(1+@1)",ROOT.RooArgList(MH,CMS_zz4l_n_sig_1_2016))
              CMS_zz4l_alpha2_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_alpha2_1_centralValue_4mu" + recobin + "2016","(1.90093+(-0.0017352)*(@0-125))",ROOT.RooArgList(MH))
              CMS_zz4l_n2_1_centralValue_4murecobin2016 = ROOT.RooFormulaVar("CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2016","CMS_zz4l_n2_1_centralValue_4mu" + recobin + "2016","(2.7194+(0)*(@0-125))",ROOT.RooArgList(MH))
              trueH = ROOT.RooDoubleCB("trueH", "DoubleCB", m, CMS_zz4l_mean_sig_1_centralValue_4murecobin2016, CMS_zz4l_sigma_sig_1_centralValue_4murecobin2016, CMS_zz4l_alpha_1_centralValue_4murecobin2016, CMS_zz4l_n_1_centralValue_4murecobin2016, CMS_zz4l_alpha2_1_centralValue_4murecobin2016, CMS_zz4l_n2_1_centralValue_4murecobin2016)

          # Wrong signal combination events
          if (channel=='4mu'):
              p1_12016 = ROOT.RooRealVar("CMS_fakeH_p1_12016","p1_12016",165.0, 145.0, 185.0)
              p3_12016 = ROOT.RooRealVar("CMS_fakeH_p3_12016","p3_12016",89.0, 84.0,94.0)
              p2_12016 = ROOT.RooFormulaVar("CMS_fakeH_p2_12016","p2_12016","0.72*@0-@1",ROOT.RooArgList(p1_12016,p3_12016))
              fakeH = ROOT.RooLandau("fakeH", "landau", m, p1_12016, p2_12016)
          if (channel=='4e'):
              p1_22016 = ROOT.RooRealVar("CMS_fakeH_p1_22016","p1_22016",165.0, 145.0, 185.0)
              p3_22016 = ROOT.RooRealVar("CMS_fakeH_p3_22016","p3_22016",89.0, 84.0,94.0)
              p2_22016 = ROOT.RooFormulaVar("CMS_fakeH_p2_22016","p2_22016","0.72*@0-@1",ROOT.RooArgList(p1_22016,p3_22016))
              fakeH = ROOT.RooLandau("fakeH", "landau", m, p1_22016, p2_22016)
          if (channel=='2e2mu'):
              p1_32016 = ROOT.RooRealVar("CMS_fakeH_p1_32016","p1_32016",165.0, 145.0, 185.0)
              p3_32016 = ROOT.RooRealVar("CMS_fakeH_p3_32016","p3_32016",89.0, 84.0,94.0)
              p2_32016 = ROOT.RooFormulaVar("CMS_fakeH_p2_32016","p2_32016","0.72*@0-@1",ROOT.RooArgList(p1_32016,p3_32016))
              fakeH = ROOT.RooLandau("fakeH", "landau", m, p1_32016, p2_32016)

    # Out of acceptance events
    # (same shape as in acceptance shape)
    out_trueH = trueH.Clone()

    # Coefficients for wrong signal combination events
    if (addfakeH):
        inc_wrongfrac_ggH=inc_wrongfrac["ggH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
        inc_wrongfrac_qqH=inc_wrongfrac["VBFH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]

    else:
        inc_wrongfrac_ggH=0.0
        inc_wrongfrac_qqH=0.0

    binfrac_wrongfrac_ggH=binfrac_wrongfrac["ggH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    binfrac_wrongfrac_qqH=binfrac_wrongfrac["VBFH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]



    #numberFake_WH = number_fake["WH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    #numberFake_ZH = number_fake["ZH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    #numberFake_ttH = number_fake["ttH125_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]

    #n_fakeH = numberFake_WH + numberFake_ZH + numberFake_ttH



    # signal shape in different recobin
    trueH_shape = {}
    fideff = {}
    fideff_var = {}
    trueH_norm = {}

    # nuisance describes the jet energy scale uncertainty
    #JES = ROOT.RooRealVar("JES","JES", 0, -5.0, 5.0)
    #if (obsName == "nJets"  or ("jet" in obsName)):
    #    lambda_JES_sig = lambdajesup[modelName+"_"+channel+"_"+obsName+"_genbin"+str(obsBin)+""+"_"+recobin]
    #    lambda_JES_sig_var = ROOT.RooRealVar("lambda_sig_"+modelName+"_"+channel+"_"+obsName+"_genbin"+str(obsBin)+""+"_"+recobin, "lambda_sig_"+modelName+"_"+channel+"_"+obsName+"_genbin"+str(obsBin)+""+"_"+recobin, lambda_JES_sig)
    #    JES_sig_rfv = ROOT.RooFormulaVar("JES_rfv_sig_"+recobin+"_"+channel,"@0*@1", ROOT.RooArgList(JES, lambda_JES_sig_var) )


    for genbin in range(nBins):
        trueH_shape[genbin] = trueH.Clone();
        trueH_shape[genbin].SetName("trueH"+channel+"Bin"+str(genbin))
        if (usecfactor): fideff[genbin] = cfactor[modelName+"_"+channel+"_"+obsName+"_genbin"+str(genbin)+"_"+recobin]
        else: fideff[genbin] = eff[modelName+"_"+channel+"_"+obsName+"_genbin"+str(genbin)+"_"+recobin]
        print "fideff[genbin]", fideff[genbin]
        print "model name is ", modelName
        fideff_var[genbin] = ROOT.RooRealVar("effBin"+str(genbin)+"_"+recobin+"_"+channel+"_"+year,"effBin"+str(genbin)+"_"+recobin+"_"+channel+"_"+year, fideff[genbin]);

        #if(not("jet" in obsName)):
        trueH_norm[genbin] = ROOT.RooFormulaVar("trueH"+channel+"Bin"+str(genbin)+"_norm","@0*@1", ROOT.RooArgList(fideff_var[genbin], lumi) );
        #else:
        #    trueH_norm[genbin] = ROOT.RooFormulaVar("trueH"+channel+"Bin"+str(genbin)+"_norm","@0*@1*(1-@2)", ROOT.RooArgList(fideff_var[genbin], lumi, JES_sig_rfv) );

    trueH_norm_final = {}
    fracBin = {}
    rBin = {}
    rBin_channel = {}
    fracSM4eBin = {}
    fracSM4muBin = {}
    K1Bin = {}
    K2Bin = {}
    SigmaBin = {}
    SigmaHBin = {}

    for genbin in range(nBins):
        if (physicalModel=="v3"):
            fid_highmass = {}
            for fState in ['4e','4mu', '2e2mu']:
                fid_highmass[fState] = 0
                fid_highmass[fState] += higgs_highmass['ggH_125.0']*higgs4l_br['125.0_'+fState]*acc['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                fid_highmass[fState] += higgs_highmass['VBF_125.0']*higgs4l_br['125.0_'+fState]*acc['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                fid_highmass[fState] += higgs_highmass['WH_125.0']*higgs4l_br['125.0_'+fState]*acc['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                fid_highmass[fState] += higgs_highmass['ZH_125.0']*higgs4l_br['125.0_'+fState]*acc['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                fid_highmass[fState] += higgs_highmass['ttH_125.0']*higgs4l_br['125.0_'+fState]*acc['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
            fid_highmass['4l'] = fid_highmass['4e'] + fid_highmass['4mu'] + fid_highmass['2e2mu']

            fracSM4eBin[str(genbin)] = ROOT.RooRealVar('fracSM4eBin'+str(genbin), 'fracSM4eBin'+str(genbin), fid_highmass['4e']/fid_highmass['4l'])
            fracSM4eBin[str(genbin)].setConstant(True)
            fracSM4muBin[str(genbin)] = ROOT.RooRealVar('fracSM4muBin'+str(genbin), 'fracSM4muBin'+str(genbin), fid_highmass['4mu']/fid_highmass['4l'])
            fracSM4muBin[str(genbin)].setConstant(True)
            K1Bin[str(genbin)] = ROOT.RooRealVar('K1Bin'+str(genbin), 'K1Bin'+str(genbin), 1.0, 0.0,  1.0/fracSM4eBin[str(genbin)].getVal())
            K2Bin[str(genbin)] = ROOT.RooRealVar('K2Bin'+str(genbin), 'K2Bin'+str(genbin), 1.0, 0.0, (1.0-fracSM4eBin[str(genbin)].getVal())/fracSM4muBin[str(genbin)].getVal())
            SigmaBin[str(genbin)] = ROOT.RooRealVar('SigmaBin'+str(genbin), 'SigmaBin'+str(genbin), fid_highmass['4l'], 0.0, 10.0)
            SigmaHBin['4e'+str(genbin)] = ROOT.RooFormulaVar("Sigma4eBin"+str(genbin),"(@0*@1*@2)", ROOT.RooArgList(SigmaBin[str(genbin)], fracSM4eBin[str(genbin)], K1Bin[str(genbin)]))
            SigmaHBin['4mu'+str(genbin)] = ROOT.RooFormulaVar("Sigma4muBin"+str(genbin),"(@0*(1.0-@1*@2)*@3*@4/(1.0-@1))", ROOT.RooArgList(SigmaBin[str(genbin)], fracSM4eBin[str(genbin)], K1Bin[str(genbin)], K2Bin[str(genbin)], fracSM4muBin[str(genbin)]))
            SigmaHBin['2e2mu'+str(genbin)] = ROOT.RooFormulaVar("Sigma2e2muBin"+str(genbin),"(@0*(1.0-@1*@2)*(1.0-@3*@4/(1.0-@1)))", ROOT.RooArgList(SigmaBin[str(genbin)], fracSM4eBin[str(genbin)], K1Bin[str(genbin)], K2Bin[str(genbin)], fracSM4muBin[str(genbin)]))
            #if ("jet" in obsName):
            #    trueH_norm_final[genbin] = ROOT.RooFormulaVar("trueH"+channel+"Bin"+str(genbin)+recobin+"_final","@0*@1*@2*(1-@3)" ,ROOT.RooArgList(SigmaHBin[channel+str(genbin)],fideff_var[genbin],lumi,JES_sig_rfv))
            #else:
            trueH_norm_final[genbin] = ROOT.RooFormulaVar("trueH"+channel+"Bin"+str(genbin)+recobin+"_final","@0*@1*@2" ,ROOT.RooArgList(SigmaHBin[channel+str(genbin)],fideff_var[genbin],lumi))
        elif (physicalModel=="v2"):
            rBin_channel[str(genbin)] = ROOT.RooRealVar("r"+channel+"Bin"+str(genbin),"r"+channel+"Bin"+str(genbin), 1.0, 0.0, 10.0)
            rBin_channel[str(genbin)].setConstant(True)
            trueH_norm_final[genbin] = ROOT.RooFormulaVar("trueH"+channel+"Bin"+str(genbin)+recobin+year+"_final","@0*@1*@2", ROOT.RooArgList(rBin_channel[str(genbin)], fideff_var[genbin],lumi))
        elif (physicalModel=="v4"):
            if channel == '2e2mu':
                rBin_channel['2e2mu'+str(genbin)] = ROOT.RooRealVar("r"+channel+"Bin"+str(genbin),"r"+channel+"Bin"+str(genbin), 1.0, 0.0, 10.0)
                trueH_norm_final[genbin] = ROOT.RooFormulaVar("trueH"+channel+"Bin"+str(genbin)+recobin+year+"_final","@0*@1*@2", ROOT.RooArgList(rBin_channel['2e2mu'+str(genbin)], fideff_var[genbin],lumi))
            else: #4e+4mu
                fid_highmass = {}
                for fState in ['4e','4mu']:
                    fid_highmass[fState] = 0
                    fid_highmass[fState] += higgs_highmasss['ggH_125.0']*higgs4l_br['125.0_'+fState]*acc['ggH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                    fid_highmass[fState] += higgs_highmasss['VBF_125.0']*higgs4l_br['125.0_'+fState]*acc['VBFH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                    fid_highmass[fState] += higgs_highmasss['WH_125.0']*higgs4l_br['125.0_'+fState]*acc['WH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                    fid_highmass[fState] += higgs_highmasss['ZH_125.0']*higgs4l_br['125.0_'+fState]*acc['ZH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                    fid_highmass[fState] += higgs_highmass['ttH_125.0']*higgs4l_br['125.0_'+fState]*acc['ttH125_'+fState+'_'+obsName+'_genbin'+str(genbin)+'_recobin'+str(genbin)]
                fid_highmass['4l'] = fid_highmass['4e'] + fid_highmass['4mu']
                fracSM4eBin[str(genbin)] = ROOT.RooRealVar('fracSM4eBin'+str(genbin), 'fracSM4eBin'+str(genbin), fid_highmass['4e']/fid_highmass['4l'])
                fracSM4eBin[str(genbin)].setConstant(True)
                fracSM4muBin[str(genbin)] = ROOT.RooRealVar('fracSM4muBin'+str(genbin), 'fracSM4muBin'+str(genbin), fid_highmass['4mu']/fid_highmass['4l'])
                fracSM4muBin[str(genbin)].setConstant(True)

                rBin_channel[str(genbin)] = ROOT.RooRealVar("r4lBin"+str(genbin),"r4lBin"+str(genbin), fid_highmass['4l'] , 0.0, 10.0)

                rBin_channel['4e'+str(genbin)] = ROOT.RooFormulaVar("r4eBin"+str(genbin),"(@0*@1)", ROOT.RooArgList(rBin_channel[str(genbin)], fracSM4eBin[str(genbin)]))
                rBin_channel['4mu'+str(genbin)] = ROOT.RooFormulaVar("r4muBin"+str(genbin),"(@0*@1)", ROOT.RooArgList(rBin_channel[str(genbin)], fracSM4muBin[str(genbin)]))
           # if ("jet" in obsName): # Even though we will not use v2 with jet variables, we keep this option in case of need
           #     trueH_norm_final[genbin] = ROOT.RooFormulaVar("trueH"+channel+"Bin"+str(genbin)+recobin+year+"_final","@0*@1*@2*(1-@3)", ROOT.RooArgList(rBin_channel[str(genbin)], fideff_var[genbin],lumi,JES_sig_rfv))
            #else:
                trueH_norm_final[genbin] = ROOT.RooFormulaVar("trueH"+channel+"Bin"+str(genbin)+recobin+year+"_final","@0*@1*@2", ROOT.RooArgList(rBin_channel[channel+str(genbin)], fideff_var[genbin],lumi))

    outin = outinratio[modelName+"_"+channel+"_"+obsName+"_genbin"+str(obsBin)+"_"+recobin]
    print "outin",obsBin,outin
    outin_var = ROOT.RooRealVar("outfracBin_"+recobin+"_"+channel+year,"outfracBin_"+recobin+"_"+channel+year, outin);
    outin_var.setConstant(True)
    out_trueH_norm_args = ROOT.RooArgList(outin_var)
    out_trueH_norm_func = "@0*("
    for i in range(nBins):
        out_trueH_norm_args.add(trueH_norm_final[i])
        out_trueH_norm_func = out_trueH_norm_func+"@"+str(i+1)+"+"
    out_trueH_norm_func = out_trueH_norm_func.replace(str(nBins)+"+",str(nBins)+")")
    out_trueH_norm = ROOT.RooFormulaVar("out_trueH_norm",out_trueH_norm_func,out_trueH_norm_args)

    frac_qqzz = fractionsBackground['qqzz_'+channel+'_'+obsName+'_'+recobin]
    frac_qqzz_var  = ROOT.RooRealVar("frac_qqzz_"+recobin+"_"+channel+"_"+year,"frac_qqzz_"+recobin+"_"+channel+"_"+year, frac_qqzz);

    frac_ggzz = fractionsBackground['ggzz_'+channel+'_'+obsName+'_'+recobin]
    frac_ggzz_var = ROOT.RooRealVar("frac_ggzz_"+recobin+"_"+channel+"_"+year,"frac_ggzz_"+recobin+"_"+channel+"_"+year, frac_ggzz);

    # frac_zjets = fractionsBackground['ZJetsCR_AllChans_'+obsName+'_'+recobin]
    # frac_zjets_var = ROOT.RooRealVar("frac_zjet_"+recobin+"_"+channel+"_"+year,"frac_zjet_"+recobin+"_"+channel+"_"+year, frac_zjets);
    frac_zjets = fractionsBackground['ZJetsCR_'+channel+'_'+obsName+'_'+recobin]
    frac_zjets_var = ROOT.RooRealVar("frac_zjet_"+recobin+"_"+channel+"_"+year,"frac_zjet_"+recobin+"_"+channel+"_"+year, frac_zjets);


    print obsBin,"frac_qqzz",frac_qqzz,"frac_ggzz",frac_ggzz,"frac_zjets",frac_zjets
    '''
    if (obsName=="nJets" or ("jet" in obsName)):
        #######
        lambda_JES_qqzz = 0.0 #lambda_qqzz_jes[modelName+"_"+channel+"_nJets_"+recobin]
        lambda_JES_qqzz_var = ROOT.RooRealVar("lambda_qqzz_"+recobin+"_"+channel,"lambda_"+recobin+"_"+channel, lambda_JES_qqzz)
        JES_qqzz_rfv = ROOT.RooFormulaVar("JES_rfv_qqzz_"+recobin+"_"+channel,"@0*@1", ROOT.RooArgList(JES, lambda_JES_qqzz_var) )
        ####
        lambda_JES_ggzz = 0.0 #lambda_ggzz_jes[modelName+"_"+channel+"_nJets_"+recobin]
        lambda_JES_ggzz_var = ROOT.RooRealVar("lambda_ggzz_"+recobin+"_"+channel,"lambda_"+recobin+"_"+channel, lambda_JES_ggzz)
        JES_ggzz_rfv = ROOT.RooFormulaVar("JES_rfv_ggzz_"+recobin+"_"+channel,"@0*@1", ROOT.RooArgList(JES, lambda_JES_ggzz_var) )
        ####
        lambda_JES_zjets = 0.0 #lambda_zjets_jes[modelName+"_"+channel+"_nJets_"+recobin]
        lambda_JES_zjets_var = ROOT.RooRealVar("lambda_zjets_"+recobin+"_"+channel,"lambda_zjets_"+recobin+"_"+channel, lambda_JES_zjets)
        JES_zjets_rfv = ROOT.RooFormulaVar("JES_rfv_zjets_"+recobin+"_"+channel,"@0*@1", ROOT.RooArgList(JES, lambda_JES_zjets_var) )
    '''

    os.chdir('../../templates/'+year+"/"+obsName+"/")

    if doubleDiff and decimal[obsName]:
        template_qqzzName = "XSBackground_qqzz_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high)+".root"
        template_ggzzName = "XSBackground_ggzz_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high)+".root"
        template_zjetsName = "XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high)+".root"
    elif doubleDiff:
        template_qqzzName = "XSBackground_qqzz_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high))+".root"
        template_ggzzName = "XSBackground_ggzz_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high))+".root"
        template_zjetsName = "XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high))+".root"
    elif decimal[obsName]:
        template_qqzzName = "XSBackground_qqzz_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+".root"
        template_ggzzName = "XSBackground_ggzz_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+".root"
        template_zjetsName = "XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+".root"
    elif not doubleDiff:
        template_qqzzName = "XSBackground_qqzz_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+".root"
        template_ggzzName = "XSBackground_ggzz_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+".root"
        template_zjetsName = "XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+".root"
    # if (not obsName=="mass4l"):
    #     template_zjetsName = "/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/templates/"+year+"/"+obsName+"/XSBackground_ZJetsCR_AllChans_"+obsName+"_"+obsBin_low+"_"+obsBin_high+".root"
    # else:
    #     template_zjetsName = "/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/templates/"+year+"/"+obsName+"/XSBackground_ZJetsCR_"+channel+"_"+obsName+"_"+obsBin_low+"_"+obsBin_high+".root"

    qqzzTempFile = ROOT.TFile(template_qqzzName,"READ")
    if decimal[obsName] and doubleDiff: qqzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high))
    elif decimal[obsName]: qqzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high))
    elif not doubleDiff: qqzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high)))
    elif doubleDiff: qqzzTemplate = qqzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high)))
    print qqzzTempFile
    print qqzzTemplate.GetName()
    print 'qqZZ bins',qqzzTemplate.GetNbinsX(),qqzzTemplate.GetBinLowEdge(1),qqzzTemplate.GetBinLowEdge(qqzzTemplate.GetNbinsX()+1)

    ggzzTempFile = ROOT.TFile(template_ggzzName,"READ")
    if decimal[obsName] and doubleDiff: ggzzTemplate = ggzzTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high))
    elif decimal[obsName]: ggzzTemplate = ggzzTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high))
    elif not doubleDiff: ggzzTemplate = ggzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high)))
    elif doubleDiff: ggzzTemplate = ggzzTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high)))
    print 'ggZZ bins',ggzzTemplate.GetNbinsX(),ggzzTemplate.GetBinLowEdge(1),ggzzTemplate.GetBinLowEdge(ggzzTemplate.GetNbinsX()+1)

    zjetsTempFile = ROOT.TFile(template_zjetsName,"READ")
    if decimal[obsName] and doubleDiff: zjetsTemplate = zjetsTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high)+"_"+str(obsBin_2nd_low)+"_"+str(obsBin_2nd_high))
    elif decimal[obsName]: zjetsTemplate = zjetsTempFile.Get("m4l_"+obsName+"_"+str(obsBin_low)+"_"+str(obsBin_high))
    elif not doubleDiff: zjetsTemplate = zjetsTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high)))
    elif doubleDiff: zjetsTemplate = zjetsTempFile.Get("m4l_"+obsName+"_"+str(trunc(obsBin_low))+"_"+str(trunc(obsBin_high))+"_"+str(trunc(obsBin_2nd_low))+"_"+str(trunc(obsBin_2nd_high)))
    print 'zjets bins',zjetsTemplate.GetNbinsX(),zjetsTemplate.GetBinLowEdge(1),zjetsTemplate.GetBinLowEdge(zjetsTemplate.GetNbinsX()+1)

    os.chdir('../../../datacard/datacard_'+year)

    binscale = 3
    qqzzTemplateNew = ROOT.TH1F("qqzzTemplateNew","qqzzTemplateNew",binscale*qqzzTemplate.GetNbinsX(),qqzzTemplate.GetBinLowEdge(1),qqzzTemplate.GetBinLowEdge(qqzzTemplate.GetNbinsX()+1))
    for i in range(1,qqzzTemplate.GetNbinsX()+1):
        for j in range(binscale):
            qqzzTemplateNew.SetBinContent((i-1)*binscale+j+1,qqzzTemplate.GetBinContent(i)/binscale)
    ggzzTemplateNew = ROOT.TH1F("ggzzTemplateNew","ggzzTemplateNew",binscale*ggzzTemplate.GetNbinsX(),ggzzTemplate.GetBinLowEdge(1),ggzzTemplate.GetBinLowEdge(ggzzTemplate.GetNbinsX()+1))
    for i in range(1,ggzzTemplate.GetNbinsX()+1):
        for j in range(binscale):
            ggzzTemplateNew.SetBinContent((i-1)*binscale+j+1,ggzzTemplate.GetBinContent(i)/binscale)
    zjetsTemplateNew = ROOT.TH1F("zjetsTemplateNew","zjetsTemplateNew",binscale*zjetsTemplate.GetNbinsX(),zjetsTemplate.GetBinLowEdge(1),zjetsTemplate.GetBinLowEdge(zjetsTemplate.GetNbinsX()+1))
    for i in range(1,zjetsTemplate.GetNbinsX()+1):
        for j in range(binscale):
            zjetsTemplateNew.SetBinContent((i-1)*binscale+j+1,zjetsTemplate.GetBinContent(i)/binscale)

    qqzzTemplateName = "qqzz_"+channel+recobin+year
    ggzzTemplateName = "ggzz_"+channel+recobin+year
    zjetsTemplateName = "zjets_"+channel+recobin+year

    qqzzTempDataHist = ROOT.RooDataHist(qqzzTemplateName,qqzzTemplateName,ROOT.RooArgList(m),qqzzTemplateNew)
    ggzzTempDataHist = ROOT.RooDataHist(ggzzTemplateName,ggzzTemplateName,ROOT.RooArgList(m),ggzzTemplateNew)
    zjetsTempDataHist = ROOT.RooDataHist(zjetsTemplateName,zjetsTemplateName,ROOT.RooArgList(m),zjetsTemplateNew)

    qqzzTemplatePdf = ROOT.RooHistPdf("qqzz","qqzz",ROOT.RooArgSet(m),qqzzTempDataHist)
    ggzzTemplatePdf = ROOT.RooHistPdf("ggzz","ggzz",ROOT.RooArgSet(m),ggzzTempDataHist)
    zjetsTemplatePdf = ROOT.RooHistPdf("zjets","zjets",ROOT.RooArgSet(m),zjetsTempDataHist)

    # bkg fractions in reco bin; implemented in terms of fractions

    #if( not (obsName=='nJets' or ("jet" in obsName) ) or (not doJES)) :
    qqzz_norm = ROOT.RooFormulaVar("bkg_qqzz_norm", "@0", ROOT.RooArgList(frac_qqzz_var) )
    ggzz_norm = ROOT.RooFormulaVar("bkg_ggzz_norm", "@0", ROOT.RooArgList(frac_ggzz_var) )
    zjets_norm = ROOT.RooFormulaVar("bkg_zjets_norm", "@0", ROOT.RooArgList(frac_zjets_var) )
    #else :
    #    qqzz_norm = ROOT.RooFormulaVar("bkg_qqzz_norm", "@0*(1-@1)", ROOT.RooArgList(frac_qqzz_var, JES_qqzz_rfv) )
    #    ggzz_norm = ROOT.RooFormulaVar("bkg_ggzz_norm", "@0*(1-@1)", ROOT.RooArgList(frac_ggzz_var, JES_ggzz_rfv) )
    #    zjets_norm = ROOT.RooFormulaVar("bkg_zjets_norm", "@0*(1-@1)", ROOT.RooArgList(frac_zjets_var, JES_zjets_rfv) )

    # Data
    # if not os.path.isfile('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/reducedTree_'+year+'.root'):
    #     os.chdir('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS')
    #     os.system('c++ -o  skim_data_tree skim_data_tree.cpp `root-config --cflags --glibs`')
    #     os.system('./skim_data_tree '+year)
    #     os.chdir('/eos/user/a/atarabin/CMSSW_10_2_13/src/HiggsAnalysis/FiducialXS/datacard_'+year)
    data_obs_file = ROOT.TFile(path['eos_path']+'Data/reducedTree_AllData_'+year+'.root')
    data_obs_tree = data_obs_file.Get('SR')

    print obsName,obsBin_low,obsBin_high
    chan = ROOT.RooRealVar("chan", "chan", 0, 3)
    # if (obsName == "nJets"): obsName = "njets_reco_pt30_eta4p7"
    if (channel=='4mu'):
        if (obsName.startswith("mass4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && chan == 1)")
        # elif (obsName.startswith("rapidity4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>105.0 && CMS_zz4l_mass<160.0 && abs("+obsName_help+")>="+obsBin_low+" && abs("+obsName_help+")<"+obsBin_high+" && chan == 1)")
        elif not doubleDiff and not obsName.startswith("mass4l"): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && chan == 1)")
        elif doubleDiff: data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,observable_2nd,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && "+obsName_2nd_help+">="+str(obsBin_2nd_low)+" && "+obsName_2nd_help+"<"+str(obsBin_2nd_high)+" && chan == 1)")
        print data_obs.numEntries()
    if (channel=='4e'):
        if (obsName.startswith("mass4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && chan == 2)")
        # elif (obsName.startswith("rapidity4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>105.0 && CMS_zz4l_mass<160.0 && abs("+obsName_help+")>="+obsBin_low+" && abs("+obsName_help+")<"+obsBin_high+" && chan == 2)")
        elif not doubleDiff and not obsName.startswith("mass4l"): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && chan == 2)")
        elif doubleDiff: data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,observable_2nd,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && "+obsName_2nd_help+">="+str(obsBin_2nd_low)+" && "+obsName_2nd_help+"<"+str(obsBin_2nd_high)+" && chan == 2)")
        print data_obs.numEntries()
    if (channel=='2e2mu'):
        if (obsName.startswith("mass4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && chan == 3)")
        # elif (obsName.startswith("rapidity4l")): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>105.0 && CMS_zz4l_mass<160.0 && abs("+obsName_help+")>="+obsBin_low+" && abs("+obsName_help+")<"+obsBin_high+" && chan == 3)")
        elif not doubleDiff and not obsName.startswith("mass4l"): data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && chan == 3)")
        elif doubleDiff: data_obs = ROOT.RooDataSet("data_obs","data_obs",data_obs_tree,ROOT.RooArgSet(m,observable,observable_2nd,chan),"(CMS_zz4l_mass>"+str(lowerBound)+" && CMS_zz4l_mass<"+str(upperBound)+" && "+obsName_help+">="+str(obsBin_low)+" && "+obsName_help+"<"+str(obsBin_high)+" && "+obsName_2nd_help+">="+str(obsBin_2nd_low)+" && "+obsName_2nd_help+"<"+str(obsBin_2nd_high)+" && chan == 3)")
        print data_obs.numEntries()
    data_obs_file.Close()

    # Create workspace
    wout = ROOT.RooWorkspace("w","w")
    if (channel=='2e2mu'):
        if (year == '2016'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_3_centralValue_2e2murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_3_2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_3_2016,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2017'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_3_centralValue_2e2murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_3_2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_3_2017,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2018'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_3_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_3_centralValue_2e2murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_3_2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_3_2018,ROOT.RooFit.RecycleConflictNodes())
    if (channel=='4e'):
        if (year == '2016'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_2_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_2_centralValue_4erecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_2_2016,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2017'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_2_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_2_centralValue_4erecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_2_2017,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2018'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_2_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_2_centralValue_4erecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_e_err_2_2018,ROOT.RooFit.RecycleConflictNodes())

    if (channel=='4mu'):
        if (year == '2016'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_1_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_1_centralValue_4murecobin2016,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_1_2016,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2017'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_1_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_1_centralValue_4murecobin2017,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_1_2017,ROOT.RooFit.RecycleConflictNodes())
        if (year == '2018'):
            getattr(wout,'import')(CMS_zz4l_mean_sig_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_sig_NoConv_1_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_sigma_sig_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_alpha2_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_n2_1_centralValue_4murecobin2018,ROOT.RooFit.RecycleConflictNodes())
            getattr(wout,'import')(CMS_zz4l_mean_m_err_1_2018,ROOT.RooFit.RecycleConflictNodes())

    for genbin in range(nBins):
        getattr(wout,'import')(trueH_shape[genbin],ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
        getattr(wout,'import')(trueH_norm[genbin],ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())

    if (not usecfactor):
        out_trueH.SetName("out_trueH")
        getattr(wout,'import')(out_trueH,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
        getattr(wout,'import')(out_trueH_norm,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())

    getattr(wout,'import')(fakeH,ROOT.RooFit.Silence())
    getattr(wout,'import')(fakeH_norm,ROOT.RooFit.Silence())

    #print "trueH norm: ",n_trueH,"fakeH norm:",n_fakeH
    qqzzTemplatePdf.SetName("bkg_qqzz")
    qqzzTemplatePdf.Print("v")
    getattr(wout,'import')(qqzzTemplatePdf,ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
    getattr(wout,'import')(qqzz_norm,ROOT.RooFit.Silence())

    ggzzTemplatePdf.SetName("bkg_ggzz")
    ggzzTemplatePdf.Print("v")
    getattr(wout,'import')(ggzzTemplatePdf,ROOT.RooFit.RecycleConflictNodes())
    getattr(wout,'import')(ggzz_norm,ROOT.RooFit.Silence())

    zjetsTemplatePdf.SetName("bkg_zjets")
    zjetsTemplatePdf.Print("v")
    getattr(wout,'import')(zjetsTemplatePdf, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
    getattr(wout,'import')(zjets_norm,ROOT.RooFit.Silence())

    ## data
    # getattr(wout,'import')(data_obs.reduce(ROOT.RooArgSet(m)),ROOT.RooFit.Silence()) #AT https://root-forum.cern.ch/t/rooworkspace-import-roofit-silence-does-not-work-when-importing-datasets/32591
    getattr(wout,'import')(data_obs.reduce(ROOT.RooArgSet(m)))

    if (addfakeH):
        if (usecfactor):
            fout = ROOT.TFile('hzz4l_'+channel+'S_13TeV_highmass_'+modelName+'_'+obsName+'_'+physicalModel+'.Databin'+str(obsBin)+'.Cfactor.root','RECREATE')
        else:
            fout = ROOT.TFile('hzz4l_'+channel+'S_13TeV_highmass_'+modelName+'_'+obsName+'_'+physicalModel+'.Databin'+str(obsBin)+'.root','RECREATE')
    else:
        if (usecfactor):
            fout = ROOT.TFile('hzz4l_'+channel+'S_13TeV_highmass_'+modelName+'_'+obsName+'_'+physicalModel+'.Databin'+str(obsBin)+'.Cfactor.NoFakeH.root','RECREATE')
        else:
            fout = ROOT.TFile('hzz4l_'+channel+'S_13TeV_highmass_'+modelName+'_'+obsName+'_'+physicalModel+'.Databin'+str(obsBin)+'.NoFakeH.root','RECREATE')

    print "write ws to fout"
    fout.WriteTObject(wout)
    fout.Close()
    os.chdir('../../fit')
    return data_obs.numEntries()

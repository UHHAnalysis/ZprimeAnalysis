#!/usr/env python

#  Script.py
#
#
#  Created by Paul Turner on 7/1/14.
#

from ROOT import *

from glob import glob

from array import array

bflavor = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_bflavor/results/ZprimePostSelectionCycle.MC.W*Jets_bflavor*root')

cflavor = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_cflavor/results/ZprimePostSelectionCycle.MC.W*Jets_cflavor*root')

lflavor = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_lflavor/results/ZprimePostSelectionCycle.MC.W*Jets_lflavor*root')

matchingup_b = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_bflavor/results/ZprimePostSelectionCycle.MC.W*matchingup*root')

matchingdown_b = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_bflavor/results/ZprimePostSelectionCycle.MC.W*matchingdown*root')

matchingup_c = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_cflavor/results/ZprimePostSelectionCycle.MC.W*matchingup*root')

matchingdown_c = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_cflavor/results/ZprimePostSelectionCycle.MC.W*matchingdown*root')

matchingup_l = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_lflavor/results/ZprimePostSelectionCycle.MC.W*matchingup*root')

matchingdown_l = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_lflavor/results/ZprimePostSelectionCycle.MC.W*matchingdown*root')

scaleup_b = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_bflavor/results/ZprimePostSelectionCycle.MC.W*scaleup*root')

scaledown_b = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_bflavor/results/ZprimePostSelectionCycle.MC.W*scaledown*root')

scaleup_c = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_cflavor/results/ZprimePostSelectionCycle.MC.W*scaleup*root')

scaledown_c = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_cflavor/results/ZprimePostSelectionCycle.MC.W*scaledown*root')

scaleup_l = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_lflavor/results/ZprimePostSelectionCycle.MC.W*scaleup*root')

scaledown_l = glob('/uscms_data/d2/drberry/Limits/CMSSW_5_3_17/src/SFrame/ZprimeAnalysis/config/ZprimePostSelection_Electron_MC_lflavor/results/ZprimePostSelectionCycle.MC.W*scaledown*root')

openfiles = {}
def getHist(name, files):
    hist = None
    for file in files:
        if file not in openfiles:
            openfiles[file] = TFile(file)
        if hist == None:
            hist = openfiles[file].Get(name).Clone()
        else:
            hist.Add(openfiles[file].Get(name))
    return hist

hists = {}

hists['bflavor'] = {}
hists['cflavor'] = {}
hists['lflavor'] = {}

hists['bflavor']['nominal'] = {}
hists['cflavor']['nominal'] = {}
hists['lflavor']['nominal'] = {}

hists['bflavor']['matchingup'] = {}
hists['cflavor']['matchingup'] = {}
hists['lflavor']['matchingup'] = {}

hists['bflavor']['matchingdown'] = {}
hists['cflavor']['matchingdown'] = {}
hists['lflavor']['matchingdown'] = {}

hists['bflavor']['scaleup'] = {}
hists['cflavor']['scaleup'] = {}
hists['lflavor']['scaleup'] = {}

hists['bflavor']['scaledown'] = {}
hists['cflavor']['scaledown'] = {}
hists['lflavor']['scaledown'] = {}

for histname in ('Kinesel','Chi2sel','TopTag','NoTopTagBTag','NoTopTagNoBTag'):
    hists['bflavor']['nominal'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', bflavor)
    hists['cflavor']['nominal'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', cflavor)
    hists['lflavor']['nominal'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', lflavor)
    hists['bflavor']['matchingup'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', matchingup_b)
    hists['cflavor']['matchingup'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', matchingup_c)
    hists['lflavor']['matchingup'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', matchingup_l)
    hists['bflavor']['matchingdown'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', matchingdown_b)
    hists['cflavor']['matchingdown'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', matchingdown_c)
    hists['lflavor']['matchingdown'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', matchingdown_l)
    hists['bflavor']['scaleup'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', scaleup_b)
    hists['cflavor']['scaleup'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', scaleup_c)
    hists['lflavor']['scaleup'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', scaleup_l)
    hists['bflavor']['scaledown'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', scaledown_b)
    hists['cflavor']['scaledown'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', scaledown_c)
    hists['lflavor']['scaledown'][histname] = getHist('Chi2_'+histname+'/M_ttbar_rec', scaledown_l)

for flavor in hists:
    print
    print
    print flavor
    for histname in ('TopTag','NoTopTagBTag','NoTopTagNoBTag'):
        print
        print histname
        nomerror = Double(0.0)
        nomint = hists[flavor]['nominal'][histname].IntegralAndError( 0, hists[flavor]['nominal'][histname].GetNbinsX() + 1, nomerror)
        matchinguperror = Double(0.0)
        matchingupint = hists[flavor]['matchingup'][histname].IntegralAndError( 0, hists[flavor]['matchingup'][histname].GetNbinsX() + 1, matchinguperror)
        matchingdownerror = Double(0.0)
        matchingdownint = hists[flavor]['matchingdown'][histname].IntegralAndError( 0, hists[flavor]['matchingdown'][histname].GetNbinsX() + 1, matchingdownerror)
        scaleuperror = Double(0.0)
        scaleupint = hists[flavor]['scaleup'][histname].IntegralAndError( 0, hists[flavor]['scaleup'][histname].GetNbinsX() + 1, scaleuperror)
        scaledownerror = Double(0.0)
        scaledownint = hists[flavor]['scaledown'][histname].IntegralAndError( 0, hists[flavor]['scaledown'][histname].GetNbinsX() + 1, scaledownerror)
        print ' Nominal:       %.2f +/- %.2f' % (nomint, nomerror)
        print ' Matching Up:   %.2f +/- %.2f' % (matchingupint, matchinguperror)
        print ' Matching Down: %.2f +/- %.2f' % (matchingdownint, matchingdownerror)
        print ' Scale Up:      %.2f +/- %.2f' % (scaleupint, scaleuperror)
        print ' Scale Down:    %.2f +/- %.2f' % (scaledownint, scaledownerror)


def findEvenBinning(nbins,hist):
    start = 0
    end = hist.GetNbinsX()+1
    print "Requested",nbins
    bingoal = hist.Integral(start,end) / nbins
    bingoal_high = bingoal
    bingoal_low = 0
    foundgoal = false
    while not foundgoal:
        start = 0
        binning = []
        binning.append(start)
        for bin in range(start+1,end+1):
            integral = hist.Integral(start,bin)
            if integral > bingoal:
                pintegral = hist.Integral(start,bin-1)
                if abs(pintegral - bingoal) < abs(integral -bingoal):
                    binning.append(bin-1)
                    start = bin
                else:
                    binning.append(bin)
                    start = bin + 1
        if len(binning) == nbins+1:
            binning[len(binning)-1] = end
        else:
            binning.append(end)
        if len(binning)-1 == nbins: foundgoal = true
        if len(binning)-1 < nbins:
            bingoal_high = bingoal
            bingoal = (bingoal_high - bingoal_low) / 2
        else:
            bingoal = 1.1 * bingoal
        print len(binning)-1,bingoal
    actualbins = []
    for bin in binning:
        actualbins.append( hist.GetBinLowEdge(bin))
    return actualbins

def createPlusMinus(name, nomhist, plushist, nplus, minushist, nminus):
    plusbinning = findEvenBinning(nplus, plushist)
    minusbinning = findEvenBinning(nminus, minushist)
    nomhist_plusbinning = nomhist.Clone()
    nomhist_plusbinning = nomhist_plusbinning.Rebin(len(plusbinning)-1,'',array('d',plusbinning))
    nomhist_minusbinning = nomhist.Clone()
    nomhist_minusbinning = nomhist_minusbinning.Rebin(len(minusbinning)-1,'',array('d',minusbinning))
    plushist_plusbinning = plushist.Clone()
    plushist_plusbinning = plushist_plusbinning.Rebin(len(plusbinning)-1,'',array('d',plusbinning))
    minushist_minusbinning = minushist.Clone()
    minushist_minusbinning = minushist_minusbinning.Rebin(len(minusbinning)-1,'',array('d',minusbinning))
    plushist_plusbinning.Divide(nomhist_plusbinning)
    minushist_minusbinning.Divide(nomhist_minusbinning)
    plus = nomhist.Clone()
    minus = nomhist.Clone()
    curbin = 1
    for bin in range(1,plushist_plusbinning.GetNbinsX()+1):
        binstart = plushist_plusbinning.GetBinLowEdge(bin)
        binend = binstart + plushist_plusbinning.GetBinWidth(bin)
        ratio = plushist_plusbinning.GetBinContent(bin)
        for bin2 in range(curbin, plus.GetNbinsX()):
            if plus.GetBinCenter(bin2) > binend:
                curbin = bin2
                break;
            if plus.GetBinCenter(bin2) > binstart and plus.GetBinCenter(bin2) < binend:
                plus.SetBinContent(bin2, plus.GetBinContent(bin2) * ratio)
            else:
                print 'Egads! You did something wrong'
    curbin = 1
    for bin in range(1,minushist_minusbinning.GetNbinsX()+1):
        binstart = minushist_minusbinning.GetBinLowEdge(bin)
        binend = binstart + minushist_minusbinning.GetBinWidth(bin)
        ratio = minushist_minusbinning.GetBinContent(bin)
        for bin2 in range(curbin, minus.GetNbinsX()):
            if minus.GetBinCenter(bin2) > binend:
                curbin = bin2
                break;
            if minus.GetBinCenter(bin2) > binstart and minus.GetBinCenter(bin2) < binend:
                minus.SetBinContent(bin2, minus.GetBinContent(bin2) * ratio)
            else:
                print 'Egads! You did something wrong'
    nameplus = name + '__plus'
    nameminus = name + '__minus'
    plus.SetName(nameplus)
    minus.SetName(nameminus)
    return plus,minus

thetanames = { 'TopTag': 'el_1top_mttbar', 'NoTopTagBTag': 'el_0top1btag_mttbar', 'NoTopTagNoBTag': 'el_0top0btag_mttbar' }

channelnames = { 'lflavor': 'wlight', 'bflavor': 'wb', 'cflavor': 'wc'}

nbins = { 'matchingup': { 'cflavor' : { 'NoTopTagNoBTag': 5 }, 'lflavor' : { 'NoTopTagNoBTag' : 15, 'NoTopTagBTag' : 15} }, 'matchingdown': { 'cflavor' : { 'NoTopTagNoBTag': 10 }, 'lflavor' : { 'NoTopTagNoBTag' : 10} }, 'scaleup' : { 'cflavor' : { 'NoTopTagNoBTag' : 25, 'NoTopTagBTag' : 35 }, 'bflavor' : { 'NoTopTagNoBTag' : 5, 'NoTopTagBTag' : 10}, 'lflavor' : { 'NoTopTagNoBTag' : 35, 'NoTopTagBTag' : 34 } }, 'scaledown' : { 'cflavor' : { 'NoTopTagNoBTag' : 15, 'NoTopTagBTag' : 5 }, 'lflavor' : { 'NoTopTagNoBTag' : 25} } }

outfile = TFile('output.root','RECREATE')
for flavor in hists:
    print
    print
    print flavor
    for histname in ('TopTag','NoTopTagBTag','NoTopTagNoBTag'):
        print
        print histname
        nplus = 1
        nminus = 1
        if flavor in nbins['scaleup']:
            if histname in nbins['scaleup'][flavor]:
                nplus = nbins['scaleup'][flavor][histname]
        if flavor in nbins['scaledown']:
            if histname in nbins['scaledown'][flavor]:
                nminus = nbins['scaledown'][flavor][histname]
        scaleplus,scaleminus = createPlusMinus(thetanames[histname]+'__'+channelnames[flavor]+'__scale', hists[flavor]['nominal'][histname], hists[flavor]['scaleup'][histname], nplus, hists[flavor]['scaledown'][histname], nminus)
        nplus = 1
        nminus = 1
        if flavor in nbins['matchingup']:
            if histname in nbins['matchingup'][flavor]:
                nplus = nbins['matchingup'][flavor][histname]
        if flavor in nbins['matchingdown']:
            if histname in nbins['matchingdown'][flavor]:
                nminus = nbins['matchingdown'][flavor][histname]
        matchingplus,matchingminus = createPlusMinus(thetanames[histname]+'__'+channelnames[flavor]+'__matching', hists[flavor]['nominal'][histname], hists[flavor]['matchingup'][histname], nplus, hists[flavor]['matchingdown'][histname], nminus)
        outfile.cd()
        scaleplus.Write()
        scaleminus.Write()
        matchingplus.Write()
        matchingminus.Write()

outfile.Close()



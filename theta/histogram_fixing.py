#! /usr/bin/env python

import sys
sys.argv.append('-b')

import ROOT
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(000000000)
ROOT.gStyle.SetOptTitle(0)

from ROOT import TCanvas, TFile, TH1, THStack, TLegend

class hinfo:
  def __init__(self, name):
    fields = name.split('__')
    self.channel = fields[0]
    self.process = fields[1]
    self.systematic = 'None'
    self.shift = 'None'
    if len(fields) > 2:
      self.systematic = fields[2]
      self.shift = fields[3]


def name(channel, process, systematic = None, shift = None):
  if not systematic:
    return '__'.join([channel, process])
  return '__'.join([channel, process, systematic, shift])


def merge(old,new):
  if not old:
    old = new.Clone()
  else:
    old.Add(new)
  return old


import math


def fixFile(filename,channelsToSym):

    file = TFile(filename)
    keys = file.GetListOfKeys()

    histos = {}
    for key in keys:
        kinfo = hinfo(key.GetName())
        if kinfo.channel not in histos: histos[kinfo.channel] = {}
        if kinfo.process not in histos[kinfo.channel]: histos[kinfo.channel][kinfo.process] = {}
        if kinfo.systematic not in histos[kinfo.channel][kinfo.process]: histos[kinfo.channel][kinfo.process][kinfo.systematic] = {}
        histos[kinfo.channel][kinfo.process][kinfo.systematic][kinfo.shift] = key.ReadObj().Clone()

    for ch in histos:
        for proc in histos[ch]:
            for sys in histos[ch][proc]:
                if ch+'__'+proc+'__'+sys in channelsToSym:
                    print 'Found channel to fix:',ch,proc,sys
                    plus = histos[ch][proc][sys]['plus']
                    minus = histos[ch][proc][sys]['minus']
                    nom = histos[ch][proc]['None']['None']
                    for bin in range(0,nom.GetNbinsX()+2):
                        diff = (plus.GetBinContent(bin) - minus.GetBinContent(bin)) / 2.0
                        avg = (plus.GetBinContent(bin) + minus.GetBinContent(bin)) / 2.0
                        change = 0.0
                        if avg !=0:
                            change = diff/avg
                        if change > 1.0: change = 1.0
                        plus.SetBinContent(bin, nom.GetBinContent(bin) * (1.0 + change))
                        minus.SetBinContent(bin, nom.GetBinContent(bin) * (1.0 - change))


    output = TFile(filename.split('.')[0]+'_fixed.root', 'RECREATE')

    for ch in sorted(histos):
        for proc in sorted(histos[ch]):
            for sys in sorted(histos[ch][proc]):
                for dir in sorted(histos[ch][proc][sys]):
                    output.cd()
                    histos[ch][proc][sys][dir].Write()

fixFile('theta_input_sixchannel_rebinned.root', ['el_0top0btag_mttbar__wlight__matching_vjets','el_0top1btag_mttbar__wlight__matching_vjets',
                                                 'el_0top0btag_mttbar__wlight__scale_vjets','el_0top1btag_mttbar__wlight__scale_vjets'])




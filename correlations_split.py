#!/usr/bin/env python
"""____________________________________________________________________
joshua.wyatt.smith@cern.ch                  

A script to calculate correlations between different variables

Input:
Ntuples with any variables.

Usage:
./correlationStudy.py

Output:
png's.
eps's
"""
import argparse

def _get_args():
  parser = argparse.ArgumentParser(description='correlations studies between variables', 
    epilog=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  return parser.parse_args()
_get_args()

import ROOT
import glob
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import pandas
import numpy
from pandas.tools.plotting import scatter_matrix

#ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 6000;")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0);

bkg_samples =[ "410501.ttbar_nonallhad_P8.p3138",
    "enugamma.p3152.v010",
    "munugamma.p3152.v010",
    "taunugamma.p3152.v010",
    "Wenu.p3317.v010",
    "Wmunu.p3317.v010",
    "Wtaunu.p3317.v010",
    "eegamma.p3152.v010",
    "mumugamma.p3152.v010",
    "tautaugamma.HM",# Remember it's this HM match
    "Zee.p3317.v010",
    "Zmumu.p3317.v010",
    "Ztautau.p3317.v010",
    "VV.p3317.v010",
    "ttV.p3317.v010",
    "ST_other.p3138.v010",
    "ST_Wt_inclusive.p3138.v010",
     "fake_CR1_data1",]

def createMatrix(correlation_list, variables, label,option):
  correlation_type=option
  if correlation_type=="hfakes":
    pretty_label="hadronic fakes"
  if correlation_type=="efakes":
    pretty_label=r"$e \to \gamma$ fakes"
  if correlation_type=="signal":
    pretty_label="signal"
  if correlation_type=="prompt_bkg":
    pretty_label=r"prompt background"
  if correlation_type=="qcd":
    pretty_label=r"Lep-fake"

  correlation_array = numpy.asarray(correlation_list)
  correlation_array_reshaped_temp = correlation_array.reshape(len(correlation_list)/len(variables),
    len(variables))
  correlation_array_reshaped=correlation_array_reshaped_temp.transpose()
  data = pandas.DataFrame(correlation_array_reshaped, columns=variables)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  cax = ax.matshow(data, vmin=-1, vmax=1)
  fig.colorbar(cax)
  ticks = numpy.arange(0,len(variables),1)
  ax.xaxis.tick_bottom()
  ax.set_xticks(ticks)
  ax.set_yticks(ticks)
  ax.set_xticklabels(variables)
  ax.set_yticklabels(variables)
  plt.setp(plt.xticks()[1], rotation=90)
  ax.tick_params(labelsize=10)
  # Hide the right and top spines
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)
  # Only show ticks on the left and bottom spines
  ax.yaxis.set_ticks_position('left')
  ax.xaxis.set_ticks_position('bottom')
  if label=="dilepton":
    plt.text(3.9,0.0,r"\textit{\textbf{ATLAS}} internal",
            fontsize=18, color='black')
    plt.text(3.9,0.70,pretty_label,
            fontsize=15, color='black')
    plt.text(3.9,1.40,str(label),
            fontsize=15, color='black')
    plt.text(3.9,2.10,r"$\sqrt{s} = 13$ TeV, 36.1 fb$^{-1}$",
            fontsize=15, color='black')
  if label=="single lepton":
    plt.text(4.9,0.3,r"\textit{\textbf{ATLAS}} internal",
            fontsize=18, color='black')
    plt.text(4.9,1.3,pretty_label,
            fontsize=15, color='black')
    plt.text(4.9,2.3,str(label),
            fontsize=15, color='black')
    plt.text(4.9,3.3,r"$\sqrt{s} = 13$ TeV, 36.1 fb$^{-1}$",
            fontsize=15, color='black')

  plt.tight_layout()
  plt.savefig("Correlations/Matrix_"+label.replace(" ","_")+"_"+option+".png")
  plt.savefig("Correlations/Matrix_"+label.replace(" ","_")+"_"+option+".eps")

def createCorrelationPlots(variables_dict, region,label, correlation_type):
  ntuples = []
  ntuples_complete = []

  prompt_selection=ROOT.TCut("!( (ph_truthOrigin[selph_index1]==23 || ph_truthOrigin[selph_index1]==24 || ph_truthOrigin[selph_index1]==25 || ph_truthOrigin[selph_index1]==26 || ph_truthOrigin[selph_index1]==27 || ph_truthOrigin[selph_index1]==28 || ph_truthOrigin[selph_index1]==29 || ph_truthOrigin[selph_index1]==30 || ph_truthOrigin[selph_index1]==31 || ph_truthOrigin[selph_index1]==32 || ph_truthOrigin[selph_index1]==33 || ph_truthOrigin[selph_index1]==34 || ph_truthOrigin[selph_index1]==35 || ph_truthOrigin[selph_index1]==42) && ph_truthType[selph_index1] == 16 ) && !( abs(ph_mc_pid[selph_index1])==11 || ( ph_mcel_dr[selph_index1]<0.05 && ph_mcel_dr[selph_index1]>=0 ) )")

  hfake_selection=ROOT.TCut("((ph_truthOrigin[selph_index1]==23 || ph_truthOrigin[selph_index1]==24 || ph_truthOrigin[selph_index1]==25 || ph_truthOrigin[selph_index1]==26 || ph_truthOrigin[selph_index1]==27 || ph_truthOrigin[selph_index1]==28 || ph_truthOrigin[selph_index1]==29 || ph_truthOrigin[selph_index1]==30 || ph_truthOrigin[selph_index1]==31 || ph_truthOrigin[selph_index1]==32 || ph_truthOrigin[selph_index1]==33 || ph_truthOrigin[selph_index1]==34 || ph_truthOrigin[selph_index1]==35 || ph_truthOrigin[selph_index1]==42) && ph_truthType[selph_index1] == 16  && !( abs(ph_mc_pid[selph_index1])==11 || ( ph_mcel_dr[selph_index1]<0.05 && ph_mcel_dr[selph_index1]>=0 ) ))")

  efake_selection=ROOT.TCut("(abs(ph_mc_pid[selph_index1])==11 || ( ph_mcel_dr[selph_index1]<0.05 && ph_mcel_dr[selph_index1]>=0 ))")
 
  QCD_selection=ROOT.TCut("(( (ejets_2015 || ejets_2016) && (ejets_2015 && (HLT_e24_lhmedium_L1EM20VH || HLT_e60_lhmedium || HLT_e120_lhloose)) || (ejets_2016 && ((HLT_e26_lhtight_nod0_ivarloose && Alt$(el_pt[0] < 61000.,0)) || ((HLT_e60_lhmedium_nod0 || HLT_e140_lhloose_nod0) && Alt$(el_pt[0] > 61000,0))))) || ( (mujets_2015 || mujets_2016) && (mujets_2015 && (HLT_mu20_iloose_L1MU15 || HLT_mu50)) || (mujets_2016 && ((HLT_mu24 && Alt$(mu_pt[0] < 51000.,0)) || (HLT_mu50 && Alt$(mu_pt[0] > 51000.,0))))))")


  weight =ROOT.TCut("weight_mc*weight_pileup*ph_SF_eff[selph_index1]*ph_SF_iso[selph_index1]*weight_leptonSF*weight_jvt*weight_bTagSF_Continuous*event_norm2 * event_lumi * ph_kfactor_overall[selph_index1]")
  # A very long and horrible cut string for each channel
  if label=="single lepton":
    cut = ROOT.TCut("(((ejets_2015 || ejets_2016) && selph_index1 >= 0 && event_ngoodphotons==1 && event_njets >= 4 && event_nbjets77 >= 1 && abs(ph_mgammalept[selph_index1] - 91188) > 5000 && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1])||((mujets_2015 || mujets_2016) && Alt$(mu_pt > 27500,0) && selph_index1 >= 0 && event_ngoodphotons==1 && event_njets >= 4 && event_nbjets77 >= 1 && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1]))")
  if label=="dilepton":
    cut = ROOT.TCut("(((ee_2015 || ee_2016) && selph_index1 >= 0  && event_ngoodphotons == 1 && event_njets >=2 && event_nbjets77 >= 1 && met_met > 30000 && (event_mll < 85000 || event_mll > 95000) && (ph_mgammaleptlept[selph_index1] < 85000 || ph_mgammaleptlept[selph_index1] > 95000) && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1])||((mumu_2015 || mumu_2016) && Alt$(mu_pt > 27500,0) && selph_index1 >=0 && event_ngoodphotons == 1 && event_njets >=2 && event_nbjets77 >= 1 && met_met > 30000 && (event_mll < 85000 || event_mll > 95000) && (ph_mgammaleptlept[selph_index1] < 85000 || ph_mgammaleptlept[selph_index1] > 95000) && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1])||((emu_2015 || emu_2016) && ((Alt$(mu_pt>27500,0) && Alt$(mu_pt>el_pt,0)) || (Alt$(el_pt>mu_pt,0))) && selph_index1>=0 && event_ngoodphotons == 1 && event_njets >=2 && event_nbjets77 >= 1 && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1]))")

  # Beautify the labels
  variables=[]
  variables_pretty=[]
  for k,v in variables_dict.iteritems():
    variables.append(k)
    variables_pretty.append(v)
  for r in region:
    path1 = "/eos/atlas/atlascerngroupdisk/phys-top/toproperties/ttgamma/v010_march18/SR1S/"+r+"/"
    path2 = "/eos/atlas/atlascerngroupdisk/phys-top/toproperties/ttgamma/v010_march18/SR1/"+r+"/"
    path3 = "/eos/atlas/atlascerngroupdisk/phys-top/toproperties/ttgamma/v010_march18/QE2_cut/"+r+"/"
    if not os.path.exists(path1):
      print "Error path1! EOS is probably at fault..."
      return
    if not os.path.exists(path2):
      print "Error path2! EOS is probably at fault..."
      return
    if not os.path.exists(path3):
      print "Error path3! EOS is probably at fault..."
      return

    if not os.path.exists("Correlations"):
      os.makedirs("Correlations")
    ntuples1 = glob.glob(path1+"*")
    ntuples2 = glob.glob(path2+"*")
    ntuples3 = glob.glob(path3+"*")
    ntuples_r = ntuples1+ntuples2+ntuples3
    ntuples.append(ntuples_r)
  for l in range(0, len(ntuples)):
    for j in ntuples[l]:
      if ".root" not in j: continue
      ntuples_complete.append(j)

  if correlation_type=="qcd":
    eChain = ROOT.TChain("nominal_Loose")
  else:
    eChain = ROOT.TChain("nominal")

  other_prompt_samples = ["ttV", "ST","enugamma", "munugamma", "taunugamma", "eegamma", "mumugamma", "tautaugamma", "VV" ]
  hfake_efake_samples = ["Wenu", "Wmunu", "Wtaunu", "Zee", "Zmumu", "Ztautau", "ttbar"]
  qcd_samples = ["fake_CR1_data15_p3315_v10","fake_CR1_data16_p3315_v10"] 

  for i in ntuples_complete:
    f = ROOT.TFile(i,"r")
    if ("fake_CR1" not in i and f.GetListOfKeys().Contains("nominal")) or ("fake_CR1" in i and f.GetListOfKeys().Contains("nominal_Loose")): 
      if correlation_type=="signal":
        if "ttgamma_noallhad.p3152.v010" in i:
          print "Adding ", i, " to signal list"
          eChain.Add(i)

      if any(x in i for x in bkg_samples):
        if any(word in i for word in hfake_efake_samples): 
            if correlation_type=="efakes" or correlation_type=="hfakes":
              print "Adding ", i, " to hfake/efake chain"
              eChain.Add(i)

        elif any(word in i for word in other_prompt_samples):
          if correlation_type=="prompt_bkg":
            eChain.Add(i)
            print "Adding ", i, " to prompt background chain"

        elif any(word in i for word in qcd_samples):
          if correlation_type=="qcd":
            eChain.Add(i)
            print "Adding ", i, " to qcd chain"
        else:
          print("Don't know what overlap removal to use!")
    else: 
      continue

  c1 = ROOT.TCanvas("c1","test",10,10,800,600)
  correlations = []
  bookkeeper=[]

  print "------- Correlations -------"
  for i in range(0,len(variables)):
    for j in range(0,len(variables)):
      var_compare = variables[i]+"_"+variables[j]

      comparison = variables[j]+":"+variables[i]
      # Speed things up and only do half
      comparison_reversed = variables[i]+":"+variables[j]
      if comparison not in bookkeeper or comparison_reversed not in bookkeeper:
        bookkeeper.append(comparison)
        bookkeeper.append(comparison_reversed)

        if correlation_type=="signal":
          overlap_removal=prompt_selection
        if correlation_type=="hfakes":
          overlap_removal=hfake_selection
        if correlation_type=="efakes":
          overlap_removal=efake_selection
        if correlation_type=="prompt_bkg":
          overlap_removal=prompt_selection
        if correlation_type=="qcd":
          overlap_removal=QCD_selection
          weight = ROOT.TCut("weights_mm_ejets[17] + weights_mm_mujets[75]")

        eChain.Draw(comparison+">>+"+var_compare,weight*cut*overlap_removal,"COL2")
        h_var_compare=ROOT.gPad.GetPrimitive(var_compare)
        correlations.append(h_var_compare.GetCorrelationFactor())
        print var_compare+" = ",h_var_compare.GetCorrelationFactor()
      else:
        correlations.append(float("NaN"))

      #Uncomment if you want ALOT of pngs. Use carefully. WIP
      #c1.Print("Correlations/"+variables[i]+"_"+variables[j]+".png")

  del eChain

  createMatrix(correlation_list=correlations, 
    variables=variables_pretty, label=label, option=correlation_type)


# createCorrelationPlots({dict of variables},[list of channels],label)
from collections import OrderedDict
var_inputs_SL=OrderedDict()
var_inputs_DL=OrderedDict()

var_inputs_SL['event_ELD_MVA_correct[selph_index1]']="ELD"
var_inputs_SL['ph_HFT_MVA[selph_index1]']="PPT"
var_inputs_SL["event_HT"]="HT"
var_inputs_SL["event_mwt"]=r"$m_{T}(W)$"
var_inputs_SL["met_met"]="MET"
var_inputs_SL["ph_mgammalept[selph_index1]"]=r"$m(\gamma,l)$"
var_inputs_SL["ph_HFT_MVA[selph_index1]"]="PPT"
var_inputs_SL["jet_pt_1st"]=r"1st jet $p_{T}$"
var_inputs_SL["jet_pt_2nd"]=r"2nd jet $p_{T}$"
var_inputs_SL["jet_pt_3rd"]=r"3rd jet $p_{T}$"
var_inputs_SL["jet_pt_4th"]=r"4th jet $p_{T}$"
var_inputs_SL["jet_pt_5th"]=r"5th jet $p_{T}$"
var_inputs_SL["jet_tagWeightBin_leading_correct"]="1st jet btag weight"
var_inputs_SL["jet_tagWeightBin_subleading_correct"]="2nd jet btag weight"
var_inputs_SL["jet_tagWeightBin_subsubleading_correct"]="3rd jet btag weight"
var_inputs_SL["event_njets"]="nr. of jets"
var_inputs_SL["event_nbjets77"]="nr. btagged jets"
#Differential variables
var_inputs_SL["ph_pt[selph_index1]"]="$p_{T}(\gamma)$"
var_inputs_SL["ph_eta[selph_index1]"]="$\eta(\gamma)$"
var_inputs_SL["ph_drlept[selph_index1]"]="$\Delta R(\gamma,l)$"

createCorrelationPlots(var_inputs_SL,
  ["ejets","mujets"], label="single lepton",correlation_type="signal")
createCorrelationPlots(var_inputs_SL,
  ["ejets","mujets"], label="single lepton",correlation_type="hfakes")
createCorrelationPlots(var_inputs_SL,
  ["ejets","mujets"], label="single lepton",correlation_type="efakes")
createCorrelationPlots(var_inputs_SL,
  ["ejets","mujets"], label="single lepton",correlation_type="prompt_bkg")
createCorrelationPlots(var_inputs_SL,
  ["ejets","mujets"], label="single lepton",correlation_type="qcd")

var_inputs_DL['event_ELD_MVA_correct[selph_index1]']="ELD"
var_inputs_DL["met_met"]="MET"
var_inputs_DL["event_mll"]=r"$m(l,l)$"
var_inputs_DL["jet_pt_1st"]=r"1st jet $p_{T}$"
var_inputs_DL["jet_pt_2nd"]=r"2nd jet $p_{T}$"
var_inputs_DL["jet_tagWeightBin_leading_correct"]="1st jet btag weight"
var_inputs_DL["jet_tagWeightBin_subleading_correct"]="2nd jet btag weight"
var_inputs_DL["event_nbjets77"]="nr. btagged jets"
#Differential variables
var_inputs_DL["ph_pt[selph_index1]"]="$p_{T}(\gamma)$"
var_inputs_DL["ph_eta[selph_index1]"]="$\eta(\gamma)$"
var_inputs_DL["ph_drlept[selph_index1]"]="$\Delta R(\gamma,l)$"
var_inputs_DL["dPhi_lep"]="$\Delta \phi(l,l)$"
var_inputs_DL["dEta_lep"]="$\Delta \eta(l,l)$"

#createCorrelationPlots(var_inputs_DL,
#  ["ee","emu","mumu"], label="dilepton", correlation_type="signal")

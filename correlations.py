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
    "ST_Wt_inclusive.p3138.v010"]

def createMatrix(correlation_list, variables, label,option):
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
    plt.text(3.9,0.70,option,
            fontsize=15, color='black')
    plt.text(3.9,1.40,str(label),
            fontsize=15, color='black')
    plt.text(3.9,2.10,r"$\sqrt{s} = 13$ TeV, 36.1 fb$^{-1}$",
            fontsize=15, color='black')
  if label=="single lepton":
    plt.text(3.9,0.0,r"\textit{\textbf{ATLAS}} internal",
            fontsize=18, color='black')
    plt.text(3.9,0.70,option,
            fontsize=15, color='black')
    plt.text(3.9,1.4,str(label),
            fontsize=15, color='black')
    plt.text(3.9,2.10,r"$\sqrt{s} = 13$ TeV, 36.1 fb$^{-1}$",
            fontsize=15, color='black')

  plt.tight_layout()
  plt.savefig("Correlations/Matrix_"+label.replace(" ","_")+"_"+option+".png")
  plt.savefig("Correlations/Matrix_"+label.replace(" ","_")+"_"+option+".eps")

def createCorrelationPlots(variables_dict, region,label):
  ntuples = []
  ntuples_complete = []

  prompt_selection="!( (ph_truthOrigin[selph_index1]==23 || ph_truthOrigin[selph_index1]==24 || ph_truthOrigin[selph_index1]==25 || ph_truthOrigin[selph_index1]==26 || ph_truthOrigin[selph_index1]==27 || ph_truthOrigin[selph_index1]==28 || ph_truthOrigin[selph_index1]==29 || ph_truthOrigin[selph_index1]==30 || ph_truthOrigin[selph_index1]==31 || ph_truthOrigin[selph_index1]==32 || ph_truthOrigin[selph_index1]==33 || ph_truthOrigin[selph_index1]==34 || ph_truthOrigin[selph_index1]==35 || ph_truthOrigin[selph_index1]==42) && ph_truthType[selph_index1] == 16 ) && !( abs(ph_mc_pid[selph_index1])==11 || ( ph_mcel_dr[selph_index1]<0.05 && ph_mcel_dr[selph_index1]>=0 ) )"

  hfake_selection="((ph_truthOrigin[selph_index1]==23 || ph_truthOrigin[selph_index1]==24 || ph_truthOrigin[selph_index1]==25 || ph_truthOrigin[selph_index1]==26 || ph_truthOrigin[selph_index1]==27 || ph_truthOrigin[selph_index1]==28 || ph_truthOrigin[selph_index1]==29 || ph_truthOrigin[selph_index1]==30 || ph_truthOrigin[selph_index1]==31 || ph_truthOrigin[selph_index1]==32 || ph_truthOrigin[selph_index1]==33 || ph_truthOrigin[selph_index1]==34 || ph_truthOrigin[selph_index1]==35 || ph_truthOrigin[selph_index1]==42) && ph_truthType[selph_index1] == 16  && !( abs(ph_mc_pid[selph_index1])==11 || ( ph_mcel_dr[selph_index1]<0.05 && ph_mcel_dr[selph_index1]>=0 ) ))"

  efake_selection="(abs(ph_mc_pid[selph_index1])==11 || ( ph_mcel_dr[selph_index1]<0.05 && ph_mcel_dr[selph_index1]>=0 ))"
 
  QCD_selection="(( (ejets_2015 || ejets_2016) && (ejets_2015 && (HLT_e24_lhmedium_L1EM20VH || HLT_e60_lhmedium || HLT_e120_lhloose)) || (ejets_2016 && ((HLT_e26_lhtight_nod0_ivarloose && Alt$(el_pt[0] < 61000.,0)) || ((HLT_e60_lhmedium_nod0 || HLT_e140_lhloose_nod0) && Alt$(el_pt[0] > 61000,0))))) || ( (mujets_2015 || mujets_2016) && (mujets_2015 && (HLT_mu20_iloose_L1MU15 || HLT_mu50)) || (mujets_2016 && ((HLT_mu24 && Alt$(mu_pt[0] < 51000.,0)) || (HLT_mu50 && Alt$(mu_pt[0] > 51000.,0))))))"

  # Beautify the labels
  variables=[]
  variables_pretty=[]
  for k,v in variables_dict.iteritems():
    variables.append(k)
    variables_pretty.append(v)
  for r in region:
    path1 = "/eos/atlas/atlascerngroupdisk/phys-top/toproperties/ttgamma/v010_march18/SR1S/"+r+"/"
    path2 = "/eos/atlas/atlascerngroupdisk/phys-top/toproperties/ttgamma/v010_march18/SR1/"+r+"/"
    if not os.path.exists(path1):
      print "Error path1! EOS is probably at fault..."
      return
    if not os.path.exists(path2):
      print "Error path2! EOS is probably at fault..."
      return
    if not os.path.exists("Correlations"):
      os.makedirs("Correlations")
    ntuples1 = glob.glob(path1+"*")
    ntuples2 = glob.glob(path2+"*")
    ntuples_r = ntuples1+ntuples2
    ntuples.append(ntuples_r)
  for l in range(0, len(ntuples)):
    for j in ntuples[l]:
      if ".root" not in j: continue
      ntuples_complete.append(j)

  signal_list = ROOT.TList()
  background_list = ROOT.TList()

  other_prompt_samples = ["ttV", "ST","enugamma", "munugamma", "taunugamma", "eegamma", "mumugamma", "tautaugamma", "VV" ]
  hfake_efake_samples = ["Wenu", "Wmunu", "Wtaunu", "Zee", "Zmumu", "Ztautau", "ttbar"]

  for i in ntuples_complete:
    f = ROOT.TFile(i,"r")
    if f.GetListOfKeys().Contains("nominal"): 
      tree = f.Get("nominal")
      ROOT.gROOT.cd()

      if "ttgamma_noallhad.p3152.v010" in i:
        print "Adding ", i, " to signal list"
        ttgamma_tree_copy = tree.CopyTree(prompt_selection) 
        ttgamma_tree_copy.SetName("nominal")
        signal_list.Add(ttgamma_tree_copy)
      if any(x in i for x in bkg_samples):
        
        if any(word in i for word in hfake_efake_samples): 
          print "Adding ", i, " to hfake/efake trees"
          hfake_tree = tree.CopyTree(hfake_selection)
          hfake_tree.SetName("nominal")
          background_list.Add(hfake_tree)
          efake_tree = tree.CopyTree(efake_selection)
          efake_tree.SetName("nominal")
          background_list.Add(efake_tree)

        elif any(word in i for word in other_prompt_samples):
          prompt_bkg_tree = tree.CopyTree(prompt_selection)
          prompt_bkg_tree.SetName("nominal")
          background_list.Add(prompt_bkg_tree)
          print "Adding ", i, " to prompt background trees"
        else:
          print("Don't know what overlap removal to use!")

    else: 
      continue

  sChain = ROOT.TTree.MergeTrees(signal_list)
  bChain = ROOT.TTree.MergeTrees(background_list)
  print "written new TTrees"

  c1 = ROOT.TCanvas("c1","test",10,10,800,600)
  correlation_bkg = []
  correlation_sig = []
  bookkeeper=[]
  weight = ROOT.TCut("weight_mc*weight_pileup*ph_SF_eff[selph_index1]*ph_SF_iso[selph_index1]*weight_leptonSF*weight_jvt*weight_bTagSF_Continuous*event_norm2 * event_lumi * ph_kfactor_overall[selph_index1]")
  # A very long and horrible cut string for each channel
  if label=="single lepton":
    cut = ROOT.TCut("(((ejets_2015 || ejets_2016) && selph_index1 >= 0 && event_ngoodphotons==1 && event_njets >= 4 && event_nbjets77 >= 1 && abs(ph_mgammalept[selph_index1] - 91188) > 5000 && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1])||((mujets_2015 || mujets_2016) && Alt$(mu_pt > 27500,0) && selph_index1 >= 0 && event_ngoodphotons==1 && event_njets >= 4 && event_nbjets77 >= 1 && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1]))")
  if label=="dilepton":
    cut = ROOT.TCut("(((ee_2015 || ee_2016) && selph_index1 >= 0  && event_ngoodphotons == 1 && event_njets >=2 && event_nbjets77 >= 1 && met_met > 30000 && (event_mll < 85000 || event_mll > 95000) && (ph_mgammaleptlept[selph_index1] < 85000 || ph_mgammaleptlept[selph_index1] > 95000) && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1])||((mumu_2015 || mumu_2016) && Alt$(mu_pt > 27500,0) && selph_index1 >=0 && event_ngoodphotons == 1 && event_njets >=2 && event_nbjets77 >= 1 && met_met > 30000 && (event_mll < 85000 || event_mll > 95000) && (ph_mgammaleptlept[selph_index1] < 85000 || ph_mgammaleptlept[selph_index1] > 95000) && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1])||((emu_2015 || emu_2016) && ((Alt$(mu_pt>27500,0) && Alt$(mu_pt>el_pt,0)) || (Alt$(el_pt>mu_pt,0))) && selph_index1>=0 && event_ngoodphotons == 1 && event_njets >=2 && event_nbjets77 >= 1 && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1]))")

  print "------- Correlations -------"
  for i in range(0,len(variables)):
    for j in range(0,len(variables)):
      sig = variables[i]+"_"+variables[j]+"_sig"
      bkg = variables[i]+"_"+variables[j]+"_bkg"

      comparison = variables[j]+":"+variables[i]
      # Speed things up and only do half
      comparison_reversed = variables[i]+":"+variables[j]
      if comparison not in bookkeeper or comparison_reversed not in bookkeeper:
        bookkeeper.append(comparison)
        bookkeeper.append(comparison_reversed)

        # sig
        sChain.Draw(comparison+">>+"+sig,weight*cut,"COL2")
        h_sig=ROOT.gPad.GetPrimitive(sig)
        correlation_sig.append(h_sig.GetCorrelationFactor())
        print sig+" = ",h_sig.GetCorrelationFactor()

        # bkg
        bChain.Draw(comparison+">>"+bkg,weight*cut,"COL2")
        h_bkg=ROOT.gPad.GetPrimitive(bkg)
        correlation_bkg.append(h_bkg.GetCorrelationFactor())
        print bkg+" = ",h_bkg.GetCorrelationFactor()
      else:
        correlation_sig.append(float("NaN"))
        correlation_bkg.append(float("NaN"))

      #Uncomment if you want ALOT of pngs. Use carefully. WIP
      #c1.Print("Correlations/"+variables[i]+"_"+variables[j]+".png")

  createMatrix(correlation_list=correlation_sig, 
    variables=variables_pretty, label=label, option="signal")
  createMatrix(correlation_list=correlation_bkg, 
    variables=variables_pretty, label=label, option="background")


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

#createCorrelationPlots(var_inputs_SL,
#  ["ejets","mujets"], label="single lepton")
#   ["mujets"], label="single lepton")

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

createCorrelationPlots(var_inputs_DL,
  #["ee","emu","mumu"], label="dilepton")
  ["ee","emu"], label="dilepton")

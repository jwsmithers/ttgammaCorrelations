#!/usr/bin/env python
"""____________________________________________________________________
joshua.wyatt.smith@cern.ch 				     			

A script to calculate correlations between different variables

Setup:
lsetup root
lsetup "lcgenv -p LCG_85swan2 x86_64-slc6-gcc49-opt numpy" 
lsetup "lcgenv -p LCG_85swan2 x86_64-slc6-gcc49-opt matplotlib" 
lsetup "lcgenv -p LCG_85swan2 x86_64-slc6-gcc49-opt pandas" 

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
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import pandas
import numpy
from pandas.tools.plotting import scatter_matrix

#ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0);

bkg_samples = [	"410501.ttbar_nonallhad_P8.p3138",
		"enugamma",
		"munugamma",
		"taunugamma",
		"Wenu",
		"Wmunu",
		"Wtaunu",
		"eegamma",
		"mumugamma",
		"tautaugamma",
		"Zee",
		"Zmumu",
		"Ztautau",
		"VV",
		"ST_other",
		"ST_Wt_inclusive"]

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
	  plt.text(2.9,0.0,r"\textit{\textbf{ATLAS}} internal",
            fontsize=18, color='black')
          plt.text(2.9,0.50,option,
            fontsize=15, color='black')
          plt.text(2.9,1.0,str(label),
            fontsize=15, color='black')
        if label=="single lepton":
          plt.text(3.9,0.3,r"\textit{\textbf{ATLAS}} internal",
            fontsize=18, color='black')
          plt.text(3.9,2.00,option,
            fontsize=15, color='black')
          plt.text(3.9,3.7,str(label),
            fontsize=15, color='black')

	plt.tight_layout()
	plt.savefig("Correlations/Matrix_"+label.replace(" ","_")+"_"+option+".png")
	plt.savefig("Correlations/Matrix_"+label.replace(" ","_")+"_"+option+".eps")

def createCorrelationPlots(variables_dict, region,label):
        ntuples = []
        ntuples_complete = []
        sChain = ROOT.TChain("nominal")
        bChain = ROOT.TChain("nominal")

	# Beautify the labels
	variables=[]
	variables_pretty=[]
	for k,v in variables_dict.iteritems():
		variables.append(k)
		variables_pretty.append(v)

        for r in region:
		path1 = "/eos/atlas/atlascerngroupdisk/phys-top/toproperties/ttgamma/v010/SR1S/"+r+"/"
		path2 = "/eos/atlas/atlascerngroupdisk/phys-top/toproperties/ttgamma/v010/SR1/"+r+"/"
		print path2
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
            ntuples_complete.append(j)
	    print j
	for i in ntuples_complete:
		if "ttgamma_nonallhadronic.p3152" in i:
                		print i
				sChain.Add(i)
		if any(x in i for x in bkg_samples):
                		print i
				bChain.Add(i)

	c1 = ROOT.TCanvas("c1","test",10,10,800,600)
	correlation_bkg = []
	correlation_sig = []
        bookkeeper=[]
        weight= ROOT.TCut("weight_mc*weight_pileup*ph_SF_eff[selph_index1]*ph_SF_iso[selph_index1]*weight_leptonSF*weight_jvt*weight_bTagSF_Continuous*event_norm * event_lumi * ph_kfactor[selph_index1]")
	cut = ROOT.TCut("(selph_index1>=0)")
        #cut = ROOT.TCut("(ejets_2015 || ejets_2016) && selph_index1 >= 0 && event_ngoodphotons==1 && event_njets >= 4 && event_nbjets77 >= 1 && abs(ph_mgammalept[selph_index1] - 91188) > 5000 && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1]")

	print "------- Correlations -------"
	for i in range(0,len(variables)):
			for j in range(0,len(variables)):
				sig = variables[i]+"_"+variables[j]+"_sig"
				bkg = variables[i]+"_"+variables[j]+"_bkg"
	
				comparison = variables[j]+":"+variables[i]
				# WIP: speed things up
				comparison_reversed = variables[i]+":"+variables[j]
				if comparison not in bookkeeper or comparison_reversed not in bookkeeper:
					bookkeeper.append(comparison)
					bookkeeper.append(comparison_reversed)
	
					# sig
					sChain.Draw(comparison+">>"+sig,weight*cut,"COL2")
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

var_inputs_SL['event_ELD_MVA[selph_index1]']="ELD"
var_inputs_SL['ph_HFT_MVA[selph_index1]']="PPT"
var_inputs_SL["event_HT"]="HT"
var_inputs_SL["event_mwt"]=r"$m_{T}(W)$"
var_inputs_SL["met_met"]="MET"
var_inputs_SL["ph_HFT_MVA[selph_index1]"]="PPT"
var_inputs_SL["jet_pt_1st"]=r"1st jet $p_{T}$"
var_inputs_SL["jet_pt_2nd"]=r"2nd jet $p_{T}$"
var_inputs_SL["jet_pt_3rd"]=r"3rd jet $p_{T}$"
var_inputs_SL["jet_pt_4th"]=r"4th jet $p_{T}$"
var_inputs_SL["jet_pt_5th"]=r"5th jet $p_{T}$"
var_inputs_SL["jet_tagWeightBin_leading"]="1st jet btag weight"
var_inputs_SL["jet_tagWeightBin_subleading"]="2nd jet btag weight"
var_inputs_SL["jet_tagWeightBin_subsubleading"]="3rd jet btag weight"
var_inputs_SL["event_njets"]="nr. of jets"
var_inputs_SL["event_nbjets77"]="nr. btagged jets"

createCorrelationPlots(var_inputs_SL,
  ["ejets"], label="single lepton")

var_inputs_DL['event_ELD_MVA[selph_index1]']="ELD"
var_inputs_DL["met_met"]="MET"
var_inputs_DL["event_mll"]=r"$m(l,l)$"
var_inputs_DL["jet_pt_1st"]=r"1st jet $p_{T}$"
var_inputs_DL["jet_pt_2nd"]=r"2nd jet $p_{T}$"
var_inputs_DL["jet_tagWeightBin_leading"]="1st jet btag weight"
var_inputs_DL["jet_tagWeightBin_subleading"]="2nd jet btag weight"
var_inputs_DL["event_nbjets77"]="nr. btagged jets"

#createCorrelationPlots(var_inputs_DL,
#  ["ee","emu","mumu"], label="dilepton")

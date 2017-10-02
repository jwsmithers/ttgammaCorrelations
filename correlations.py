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
# Correction Matrix Plot
import matplotlib.pyplot as plt
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
	correlation_array_reshaped = correlation_array.reshape(len(correlation_list)/len(variables),
		len(variables))
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
	plt.title(option+" "+str(label))
	plt.tight_layout()
	plt.savefig("Correlations/Matrix_"+label.replace(" ","_")+"_"+option+".png")

def createCorrelationPlots(variables, region,label):
        ntuples = []
        ntuples_complete = []
        sChain = ROOT.TChain("nominal")
        bChain = ROOT.TChain("nominal")
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
        weight= ROOT.TCut("weight_mc*weight_pileup*ph_SF_eff[selph_index1]*ph_SF_iso[selph_index1]*weight_leptonSF*weight_jvt*weight_bTagSF_Continuous*event_norm * event_lumi * ph_kfactor[selph_index1]")
        cut = ROOT.TCut("(ejets_2015 || ejets_2016) && selph_index1 >= 0 && event_ngoodphotons==1 && event_njets >= 4 && event_nbjets77 >= 1 && abs(ph_mgammalept[selph_index1] - 91188) > 5000 && ph_drlept[selph_index1] > 1.0 && ph_isoFCT[selph_index1]")

	print "------- Correlations -------"
	for i in range(0,len(variables)):
		for j in range(0,len(variables)):
			sig = variables[i]+"_"+variables[j]+"_sig"
			bkg = variables[i]+"_"+variables[j]+"_bkg"

			comparison = variables[i]+":"+variables[j]
			# WIP: speed things up
			#comparison_reversed = variables[j]+":"+variables[i]
                        #if comparison == comparison_reversed:
			#	correlation_sig.append(1)
			#	correlation_bkg.append(1)
			#	continue

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

			#Uncomment if you want ALOT of pngs. Use carefully. WIP
			#c1.Print("Correlations/"+variables[i]+"_"+variables[j]+".png")

	createMatrix(correlation_list=correlation_sig, 
		variables=variables, label=label, option="signal")
	createMatrix(correlation_list=correlation_bkg, 
		variables=variables, label=label, option="background")


# createCorrelationPlots([list of variables],[list of channels],label)
createCorrelationPlots([
  'event_ELD_MVA', 
  'ph_pt', 
  'event_njets'], ["ejets"], label="single lepton")


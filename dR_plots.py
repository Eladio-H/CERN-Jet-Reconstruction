import ROOT
import math
import numpy as np
import sys

def getHist(inFile, histName):
    hist = inFile.Get(histName)
    if not hist:
        print("Unable to retrieve a required histogram, please make sure that the input file was run up to the required step")
        print("The missing histogram is named:", histName)
        return None
    return hist

# LOAD HISTOGRAMS FROM INPUT (ANALYSIS) ROOT FILE
inFileName  = "truthDR_analysis.root"
outFileName = "truthDR_analysis.pdf"
inFile = ROOT.TFile.Open(inFileName,"READ")

hist_deltaR04 = getHist(inFile, "deltaR04")
hist_all_pt = getHist(inFile, "all_pt")
hist_dR_pt = getHist(inFile, "dR_pt")

hist_double_any = getHist(inFile, "double_any")
hist_all_any = getHist(inFile, "all_any")

hist_double = getHist(inFile, "double")
hist_all = getHist(inFile, "all")

hist_double_1 = getHist(inFile, "double_1")
hist_all_1 = getHist(inFile, "all_1")

hist_double_any_t = getHist(inFile, "double_any_t")
hist_all_any_t = getHist(inFile, "all_any_t")

hist_double_t = getHist(inFile, "double_t")
hist_all_t = getHist(inFile, "all_t")

hist_double_1_t = getHist(inFile, "double_1_t")
hist_all_1_t = getHist(inFile, "all_1_t")

hist_num_deltaR = getHist(inFile, "num_deltaR")
hist_all_fCPV = getHist(inFile, "all_fCPV")


# START PLOTTING
canvas = ROOT.TCanvas("canvas","canvas",1080,650)
canvas.cd()
canvas.Print(outFileName+"[")

# PLOT PROBABILITY OF THERE BEING NO JET WITHIN dR < 0.4 OF A JET WITH fCPV < 1
efficiency = ROOT.TEfficiency(hist_num_deltaR, hist_all_fCPV)
efficiency.SetTitle("Frequency of there being no other jets within dR < 0.4 for a reco jet of fCPV < 1; reco pT; Probability")
canvas.SetLogx(False)
canvas.SetLogy(False)
efficiency.SetLineColor(ROOT.kBlack)
efficiency.Draw("AP")
canvas.Print(outFileName)

# PLOT PROBABILITY OF TRUTH JET OVERLAPPING AS A FUNCTION OF pT
efficiency = ROOT.TEfficiency(hist_dR_pt, hist_all_pt)
efficiency.SetTitle("Probability of a truth jet overlapping as a function of pT; truth pT; Probability")
canvas.SetLogx(False) 
canvas.SetLogy(False)
efficiency.SetLineColor(ROOT.kBlue)
efficiency.Draw("AP")
canvas.Print(outFileName)

# SAVE PLOTS TO PDF
canvas.Print(outFileName+"]")

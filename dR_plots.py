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
inFileName  = "dR_analysis.root"
outFileName = "dR_analysis.pdf"
inFile = ROOT.TFile.Open(inFileName,"READ")

hist_all_pt = getHist(inFile, "all_pt")
hist_dR_pt = getHist(inFile, "dR_pt")

hist_double = getHist(inFile, "double")
hist_double_1 = getHist(inFile, "double_1")

hist_num_deltaR = getHist(inFile, "num_deltaR")
hist_all_fCPV = getHist(inFile, "all_fCPV")

hist_num_deltaR_1 = getHist(inFile, "num_deltaR_1")
hist_all_fCPV_1 = getHist(inFile, "all_fCPV_1")

hist_gpv = getHist(inFile, "gpv")
hist_pileup = getHist(inFile, "pileup")

hist_gpv_1 = getHist(inFile, "gpv_1")
hist_pileup_1 = getHist(inFile, "pileup_1")

hist_all_filt = getHist(inFile, "all_fCPV_filt")
hist_all_filt_1 = getHist(inFile, "all_fCPV_filt_1")

hist_double_p = getHist(inFile, "double_p")
hist_double_1_p = getHist(inFile, "double_1_p")

hist_all_incl = hist_all_fCPV.Clone()
hist_all_incl.Add(hist_all_fCPV_1)


# PLOTTING FUNCTION
def plot(hist, fit, title, colour, d2, teff, xt, yt, description):
    topPad = ROOT.TPad("topPad", "Top Pad", 0, 0.1, 1, 1.0)  # Upper 90% for the plot
    bottomPad = ROOT.TPad("bottomPad", "Bottom Pad", 0, 0.0, 1, 0.1)  # Lower 10% for the text
    topPad.SetBottomMargin(0.15)
    topPad.Draw()
    topPad.cd()

    canvas.SetLogx(False)
    canvas.SetLogy(False)
    if teff:
        hist.SetTitle(title)
        hist.SetLineColor(colour)
        hist.Draw("AP")

    elif d2:
        hist.SetTitle(title)
        hist.SetStats(0)
        hist.GetXaxis().SetTitle(xt)
        hist.GetYaxis().SetTitle(yt)
        hist.Draw("colz")

    else:
        hist.SetTitle(title)
        hist.SetStats(0)
        hist.GetXaxis().SetTitle(xt)
        hist.GetYaxis().SetTitle(yt)
        hist.SetLineColor(colour)
        hist.SetLineWidth(2)
        hist.Draw()

        if fit:
            fit_gaus(hist)


    # Write description
    canvas.cd()  
    bottomPad.SetTopMargin(0.0)  # No top margin for bottom pad
    bottomPad.Draw()  # Draw the bottom pad on the canvas
    bottomPad.cd()  # Make the bottom pad current

    # Add text in the bottom pad using TLatex
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.3)  # Larger text size for clarity
    latex.SetTextAlign(22)  # Centered text
    latex.DrawLatex(0.5, 0.85, description)

    canvas.cd()

    canvas.Print(outFileName)


# START PLOTTING
canvas = ROOT.TCanvas("canvas","canvas",1080,650)
canvas.cd()
canvas.Print(outFileName+"[")

# PLOT PROBABILITY OF THERE BEING NO JET WITHIN dR < 0.4 OF A JET WITH fCPV < 1
efficiency = ROOT.TEfficiency(hist_num_deltaR, hist_all_fCPV)
description = "Filled with fraction of: no. of reco jets within dR < 0.4 of neighbour / total no. of jets of fCPV < 1; for each reco pT bin"
title = "Frequency of there being no other jets within dR < 0.4 for a reco jet of fCPV < 1; reco pT; Probability"
plot(efficiency, False, title, ROOT.kBlack, False, True, "", "", description)

# PLOT PROBABILITY OF THERE BEING NO JET WITHIN dR < 0.4 OF A JET WITH fCPV = 1
efficiency = ROOT.TEfficiency(hist_num_deltaR_1, hist_all_fCPV_1)
description = "Filled with fraction of: no. of reco jets within dR < 0.4 of neighbour / total no. of jets of fCPV < 1; for each reco pT bin"
title = "Frequency of there being no other jets within dR < 0.4 for a reco jet of fCPV = 1; reco pT; Probability"
plot(efficiency, False, title, ROOT.kBlack, False, True, "", "", description)


# PLOT PROBABILITY OF A JET BEING GPV-MATCHED
hist_incl = hist_gpv.Clone()
hist_incl.Add(hist_gpv_1)

efficiency = ROOT.TEfficiency(hist_incl, hist_all_incl)
description = "Filled with fraction of: no. of reco jets matching a GPV jet / total no. of jets of any fCPV; for each reco pT bin"
title = "Probability of a reco jet of any fCPV match a GPV-originated jet; reco pT; Probability"
plot(efficiency, False, title, ROOT.kRed, False, True, "", "", description)

efficiency = ROOT.TEfficiency(hist_gpv, hist_all_fCPV)
description = "Filled with fraction of: no. of reco jets matching a GPV jet / total no. of jets of fCPV < 1; for each reco pT bin"
title = "Probability of a reco jet of fCPV < 1 match a GPV-originated jet; reco pT; Probability"
plot(efficiency, False, title, ROOT.kRed, False, True, "", "", description)

efficiency = ROOT.TEfficiency(hist_gpv_1, hist_all_fCPV_1)
description = "Filled with fraction of: no. of reco jets matching a GPV jet / total no. of jets of fCPV = 1; for each reco pT bin"
title = "Probability of a reco jet of fCPV = 1 match a GPV-originated jet; reco pT; Probability"
plot(efficiency, False, title, ROOT.kRed, False, True, "", "", description)

# PLOT PROBABILITY OF A JET BEING PILEUP-MATCHED
hist_incl = hist_pileup.Clone()
hist_incl.Add(hist_pileup_1)
efficiency = ROOT.TEfficiency(hist_incl, hist_all_incl)
description = "Filled with fraction of: no. of reco jets matching a pileup-originated jet / total no. of jets of any fCPV; for each reco pT bin"
title = "Probability of a reco jet of any fCPV match a pileup-originated jet; reco pT; Probability"
plot(efficiency, False, title, ROOT.kBlue, False, True, "", "", description)

efficiency = ROOT.TEfficiency(hist_pileup, hist_all_fCPV)
description = "Filled with fraction of: no. of reco jets matching a pileup-originated jet / total no. of jets of fCPV < 1; for each reco pT bin"
title = "Probability of a reco jet of fCPV < 1 match a pileup-originated jet; reco pT; Probability"
plot(efficiency, False, title, ROOT.kBlue, False, True, "", "", description)

efficiency = ROOT.TEfficiency(hist_pileup_1, hist_all_fCPV_1)
description = "Filled with fraction of: no. of reco jets matching a pileup-originated jet / total no. of jets of fCPV = 1; for each reco pT bin"
title = "Probability of a reco jet of fCPV = 1 match a pileup-originated jet; reco pT; Probability"
plot(efficiency, False, title, ROOT.kBlue, False, True, "", "", description)


# PLOT PROBABILITY OF A RECO JET BEING DOUBLE-MATCHED WITH fCPV < 1
efficiency = ROOT.TEfficiency(hist_double, hist_all_fCPV)
description = "Filled with fraction of: no. of double-matched reco jets / total no. of reco jets of fCPV < 1; for each reco pT bin"
title = "Probability of a reco jet of fCPV < 1 to be double truth-matched; reco pT; Probability"
plot(efficiency, False, title, ROOT.kBlack, False, True, "", "", description)

# PLOT PROBABILITY OF A RECO JET BEING DOUBLE-MATCHED WITH fCPV = 1
efficiency = ROOT.TEfficiency(hist_double_1, hist_all_fCPV_1)
description = "Filled with fraction of: no. of double-matched reco jets / total no. of reco jets of fCPV = 1; for each reco pT bin"
title = "Probability of a reco jet of fCPV = 1 to be double truth-matched; reco pT; Probability"
plot(efficiency, False, title, ROOT.kBlack, False, True, "", "", description)

# PLOT PROBABILITY OF GPV MATCHED JET BEING DOUBLE MATCHED
efficiency = ROOT.TEfficiency(hist_double, hist_gpv)
description = "Filled with fraction of: no. of double-matched reco jets / total no. of GPV-matched reco jets of fCPV < 1; for each reco pT bin"
title = "Probability of GPV-matched reco jet of fCPV < 1 to be double matched; reco pT; Probability"
plot(efficiency, False, title, ROOT.kRed, False, True, "", "", description)

efficiency = ROOT.TEfficiency(hist_double_1, hist_gpv_1)
description = "Filled with fraction of: no. of double-matched reco jets / total no. of GPV-matched reco jets of fCPV = 1; for each reco pT bin"
title = "Probability of GPV-matched reco jet of fCPV = 1 to be double matched; reco pT; Probability"
plot(efficiency, False, title, ROOT.kRed, False, True, "", "", description)


# PLOT PROBABILITY OF PILEUP MATCHED JET BEING DOUBLE MATCHED
efficiency = ROOT.TEfficiency(hist_double_p, hist_pileup)
description = "Filled with fraction of: no. of double-matched reco jets / total no. of pileup-matched reco jets of fCPV < 1; for each reco pT bin"
title = "Probability of pileup-matched reco jet of fCPV < 1 to be double matched; reco pT; Probability"
plot(efficiency, False, title, ROOT.kBlue, False, True, "", "", description)

efficiency = ROOT.TEfficiency(hist_double_1_p, hist_pileup_1)
description = "Filled with fraction of: no. of double-matched reco jets / total no. of pileup-matched reco jets of fCPV = 1; for each reco pT bin"
title = "Probability of pileup-matched reco jet of fCPV = 1 to be double matched; reco pT; Probability"
plot(efficiency, False, title, ROOT.kBlue, False, True, "", "", description)


# PLOT PROBABILITY OF TRUTH JET OVERLAPPING AS A FUNCTION OF pT
efficiency = ROOT.TEfficiency(hist_dR_pt, hist_all_pt)
description = "Filled with fraction of: no. of truth jets that overlap / total no. of truth jets; for each truth pT bin"
title = "Probability of a truth jet overlapping as a function of pT; truth pT; Probability"
plot(efficiency, False, title, ROOT.kViolet, False, True, "", "", description)


# SAVE PLOTS TO PDF
canvas.Print(outFileName+"]")

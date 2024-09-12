import ROOT
import math
import numpy as np
import sys

def getHist(inFile, histName):
    hist = inFile.Get(histName)
    if not hist:
        print("Unable to retrieve a required histogram, please make sure that the input file was run up to the required step")
        print("The missing histogram is named:",histName)
        return None
    return hist


# LOAD HISTOGRAMS FROM INPUT (ANALYSIS) ROOT FILE
inFileName  = "ntuple_analysis.root"
outFileName = "ntuple_analysis.pdf"

inFile = ROOT.TFile.Open(inFileName,"READ")
hist_cFPV_unf = getHist(inFile, "fCPV_unfiltered")
hist_cFPV_truth = getHist(inFile, "fCPV_truth")
hist_cFPV_pileup = getHist(inFile, "fCPV_pileup")
hist_cFPV_truth_pileup = getHist(inFile, "fCPV_truth_pileup")
hist_fCPV_truth_pileup_or = getHist(inFile, "fCPV_truth_pileup_or")

hist_pt_comparison_fCPV_any = getHist(inFile, "pt_comparison_any")
hist_pt_comparison_fCPV = getHist(inFile, "pt_comparison")
hist_pt_comparison_fCPV_1 = getHist(inFile, "pt_comparison_1")

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

# In a loop, load the all pT response histograms for all the ranges into 'hists_pt_responses'
hists_pt_responses = []
ranges = [[20.e3, 25.e3], [25.e3, 30.e3], [30.e3, 40.e3]]
for i in ranges:
    base = str(i[0]) + "-" + str(i[1])
    for_range = []

    for_range.append(getHist(inFile, "pt_response_any_" + base))
    for_range.append(getHist(inFile, "pt_response_" + base))
    for_range.append(getHist(inFile, "pt_response_1_" + base))

    for_range.append(getHist(inFile, "pt_response_any_p_" + base))
    for_range.append(getHist(inFile, "pt_response_p_" + base))
    for_range.append(getHist(inFile, "pt_response_1_p_" + base))

    for_range.append(getHist(inFile, "pt_response_any_both_" + base))
    for_range.append(getHist(inFile, "pt_response_both_" + base))
    for_range.append(getHist(inFile, "pt_response_1_both_" + base))

    hists_pt_responses.append(for_range)


# Function to fit a gaussian curve on a given plot (certain histogram)
def fit_gaus(hist):
    gf = ROOT.TF1("gf", "gaus")
    maxX = hist.GetBinCenter(hist.GetMaximumBin())
    RMS = hist.GetRMS()
    gf.SetParameters(hist.GetMaximum(), maxX, RMS, 10.0)
    hist.Fit(gf, "E", "", maxX-RMS, maxX+RMS)
    hist.Draw("pe")

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.03)

    chi2 = gf.GetChisquare()
    ndof = gf.GetNDF()
    mean = gf.GetParameter(1)
    width = gf.GetParameter(2) # or 0/2
    entries = hist.GetEntries()

    x_position = 0.64  # New x position (closer to the left side)
    y_position_start = 0.80  # New y position start (higher up)

    latex.DrawText(x_position, y_position_start, "Number of entries = " + str(entries))
    latex.DrawText(x_position, y_position_start - 0.05, "Mean = %.3f"%(mean))
    latex.DrawText(x_position, y_position_start - 0.10, "Width = %.3f"%(width))
    latex.DrawText(x_position, y_position_start - 0.15, "chi2/ndof = %.1f/%d = %.1f"%(chi2, ndof, chi2/ndof))


# START PLOTTING
canvas = ROOT.TCanvas("canvas","canvas",1300,800)
canvas.cd()
canvas.Print(outFileName+"[")
text_size = 0.05


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


# PLOT FCPV COMPARISON HISTOGRAM
canvas.SetLogx(False)
canvas.SetLogy(False)

topPad = ROOT.TPad("topPad", "Top Pad", 0, 0.1, 1, 1.0)  # Upper 90% for the plot
bottomPad = ROOT.TPad("bottomPad", "Bottom Pad", 0, 0.0, 1, 0.1)  # Lower 10% for the text
topPad.SetBottomMargin(0.15)
topPad.Draw()
topPad.cd()

hist_cFPV_unf.SetLineColor(ROOT.kRed)
hist_cFPV_unf.SetLineWidth(2)
hist_cFPV_unf.GetXaxis().SetTitle("fCPV of reco jets")
hist_cFPV_unf.GetYaxis().SetTitle("Fraction of events")
hist_cFPV_unf.GetXaxis().SetLabelSize(text_size)
hist_cFPV_unf.SetStats(0)
hist_cFPV_unf.Scale(1./hist_cFPV_unf.Integral())
hist_cFPV_unf.SetMaximum(1)
hist_cFPV_unf.Draw()
hist_cFPV_truth.SetLineColor(ROOT.kBlue)
hist_cFPV_truth.SetLineWidth(2)
hist_cFPV_truth.GetXaxis().SetTitle("fCPV of truth-matched reco jets")
hist_cFPV_truth.GetYaxis().SetTitle("Scaled number of events")
hist_cFPV_truth.GetXaxis().SetLabelSize(text_size)
hist_cFPV_truth.SetStats(0)
hist_cFPV_truth.Scale(1./hist_cFPV_truth.Integral())
hist_cFPV_truth.Draw("same")
hist_cFPV_pileup.SetLineColor(ROOT.kBlack)
hist_cFPV_pileup.SetLineWidth(2)
hist_cFPV_pileup.GetXaxis().SetTitle("fCPV of truth-pileup-matched reco jets")
hist_cFPV_pileup.GetYaxis().SetTitle("Scaled number of events")
hist_cFPV_pileup.GetXaxis().SetLabelSize(text_size)
hist_cFPV_pileup.SetStats(0)
hist_cFPV_pileup.Scale(1./hist_cFPV_pileup.Integral())
hist_cFPV_pileup.Draw("same")
hist_cFPV_truth_pileup.SetLineColor(ROOT.kGreen)
hist_cFPV_truth_pileup.SetLineWidth(2)
hist_cFPV_truth_pileup.GetXaxis().SetTitle("fCPV of truth- and pileup-matched reco jets")
hist_cFPV_truth_pileup.GetYaxis().SetTitle("Scaled number of events")
hist_cFPV_truth_pileup.GetXaxis().SetLabelSize(text_size)
hist_cFPV_truth_pileup.SetStats(0)
hist_cFPV_truth_pileup.Scale(1./hist_cFPV_truth_pileup.Integral())
hist_cFPV_truth_pileup.Draw("same")
hist_fCPV_truth_pileup_or.SetLineColor(ROOT.kOrange)
hist_fCPV_truth_pileup_or.SetLineWidth(2)
hist_fCPV_truth_pileup_or.GetXaxis().SetTitle("fCPV of truth- and pileup-matched reco jets")
hist_fCPV_truth_pileup_or.GetYaxis().SetTitle("Fraction of events")
hist_fCPV_truth_pileup_or.GetXaxis().SetLabelSize(text_size)
hist_fCPV_truth_pileup_or.SetStats(0)
hist_fCPV_truth_pileup_or.Scale(1./hist_fCPV_truth_pileup_or.Integral())
hist_fCPV_truth_pileup_or.Draw("same")
legend = ROOT.TLegend(0.15, 0.75, 0.45, 0.85)
legend.AddEntry(hist_cFPV_unf)
legend.AddEntry(hist_cFPV_truth)
legend.AddEntry(hist_cFPV_pileup)
legend.AddEntry(hist_cFPV_truth_pileup)
legend.AddEntry(hist_fCPV_truth_pileup_or)
legend.SetTextSize(0.025)
legend.SetBorderSize(0)
legend.Draw("same")

canvas.cd()  
bottomPad.SetTopMargin(0.0)  # No top margin for bottom pad
bottomPad.Draw()  # Draw the bottom pad on the canvas
bottomPad.cd()  # Make the bottom pad current

# Add text in the bottom pad using TLatex
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.3)  # Larger text size for clarity
latex.SetTextAlign(22)  # Centered text
latex.DrawLatex(0.5, 0.85, "Filled with fCPV of reco jets, categorised by to what extent they are matched")

canvas.cd()

canvas.Print(outFileName)

# PLOT LEADING JET VS SUBLEADING JET FOR DOUBLE MATCHED JETS

# Any fCPV
description = "For double matched (GPV / Pileup) reco jets of any fCPV: Filled with the leading vs subleading matched jet pT"
title = "Correlation between Leading Jet p_{T} and Subleading Jet p_{T} for any fCPV"
xt = "p{T} of Leading Jet"
yt = "p{T} of Subleading Jet"
plot(hist_pt_comparison_fCPV_any, False, title, ROOT.kBlack, True, False, xt, yt, description)

# fCPV = 1
description = "For double matched (GPV / Pileup) reco jets of fCPV = 1: Filled with the leading vs subleading matched jet pT"
title = "Correlation between Leading Jet p_{T} and Subleading Jet p_{T} for fCPV = 1"
xt = "p{T} of Leading Jet"
yt = "p{T} of Subleading Jet"
plot(hist_pt_comparison_fCPV_any, False, title, ROOT.kBlack, True, False, xt, yt, description)

# fCPV < 1
description = "For double matched (GPV / Pileup) reco jets of fCPV < 1: Filled with the leading vs subleading matched jet pT"
title = "Correlation between Leading Jet p_{T} and Subleading Jet p_{T} for fCPV < 1"
xt = "p{T} of Leading Jet"
yt = "p{T} of Subleading Jet"
plot(hist_pt_comparison_fCPV_any, False, title, ROOT.kBlack, True, False, xt, yt, description)


# PLOT pT RESPONSE PLOTS
for j in range(len(ranges)):
    group = hists_pt_responses[j]
    range_str = "Range: " + str(ranges[j][0]) + "-" + str(ranges[j][1]) + " MeV"

    either_or = group[0].Clone("either_or")
    either_or.Add(group[3])

    # ANY fCPV

    # only truth-matched
    group[0].Scale(1./group[0].Integral())
    title = "pT response: p{T} of reco jets of any fCPV / p{T} of matching truth jet - " + range_str
    description = "Filled with the pT response of only GPV-matched reco jets of any fCPV - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(group[0], True, title, ROOT.kBlack, False, False, xt, yt, description)

    # only pileup-matched
    group[3].Scale(1./group[3].Integral())
    title = "pT response: p{T} of reco jets of any fCPV / p{T} of matching inTimeTruthJet - " + range_str
    description = "Filled with the pT response of only pileup-matched reco jets of any fCPV - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(group[3], True, title, ROOT.kBlack, False, False, xt, yt, description)

    # matched
    canvas.SetLogx(False)
    either_or.Scale(1./either_or.Integral())
    title = "pT response: p{T} of reco jets of any fCPV / p{T} of matching truth jet OR inTimeTruthJet - " + range_str
    description = "Filled with the pT response of single-matched (either GPV or pileup) reco jets of any fCPV - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(either_or, True, title, ROOT.kBlack, False, False, xt, yt, description)

    # double matched
    group[6].Scale(1./group[6].Integral())
    title = "pT response: p{T} of reco jets of any fCPV / p{T} of  BOTH matching truth jets - " + range_str
    description = "Filled with the pT response of double-matched reco jets of any fCPV - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(group[6], True, title, ROOT.kBlack, False, False, xt, yt, description)

    # fCPV = 1
    either_or = group[2].Clone("either_or")
    either_or.Add(group[5])

    # only truth-matched
    group[2].Scale(1./group[2].Integral())
    title = "pT response: p{T} of reco jets of fCPV = 1 / p{T} of matching truth jet - " + range_str
    description = "Filled with the pT response of only GPV-matched reco jets of fCPV = 1 - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(group[2], True, title, ROOT.kBlack, False, False, xt, yt, description)

    # only pileup-matched
    group[5].Scale(1./group[5].Integral())
    title = "pT response: p{T} of reco jets of fCPV = 1 / p{T} of matching inTimeTruthJet - " + range_str
    description = "Filled with the pT response of only pileup-matched reco jets of fCPV = 1 - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(group[5], True, title, ROOT.kBlack, False, False, xt, yt, description)

    # matched
    canvas.SetLogx(False)
    either_or.Scale(1./either_or.Integral())
    title = "pT response: p{T} of reco jets of fCPV = 1 / p{T} of matching truth jet OR inTimeTruthJet - " + range_str
    description = "Filled with the pT response of single-matched (either GPV or pileup) reco jets of fCPV = 1 - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(either_or, True, title, ROOT.kBlack, False, False, xt, yt, description)

    # double matched
    group[8].Scale(1./group[8].Integral())
    title = "pT response: p{T} of reco jets of fCPV = 1 / p{T} of  BOTH matching truth jets - " + range_str
    description = "Filled with the pT response of double-matched reco jets of fCPV = 1 - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(group[8], True, title, ROOT.kBlack, False, False, xt, yt, description)

    # fCPV < 1
    either_or = group[1].Clone("either_or")
    either_or.Add(group[4])

    # only truth-matched
    group[1].Scale(1./group[1].Integral())
    title = "pT response: p{T} of reco jets of fCPV < 1 / p{T} of matching truth jet - " + range_str
    description = "Filled with the pT response of only GPV-matched reco jets of fCPV < 1 - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(group[1], True, title, ROOT.kBlack, False, False, xt, yt, description)

    # only pileup-matched
    group[4].Scale(1./group[4].Integral())
    title = "pT response: p{T} of reco jets of fCPV < 1 / p{T} of matching inTimeTruthJet - " + range_str
    description = "Filled with the pT response of only pileup-matched reco jets of fCPV < 1 - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(group[4], True, title, ROOT.kBlack, False, False, xt, yt, description)

    # matched
    canvas.SetLogx(False)
    either_or.Scale(1./either_or.Integral())
    title = "pT response: p{T} of reco jets of fCPV < 1 / p{T} of matching truth jet OR inTimeTruthJet - " + range_str
    description = "Filled with the pT response of single-matched (either GPV or pileup) reco jets of fCPV < 1 - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(either_or, True, title, ROOT.kBlack, False, False, xt, yt, description)

    # double matched
    group[7].Scale(1./group[7].Integral())
    title = "pT response: p{T} of reco jets of fCPV < 1 / p{T} of  BOTH matching truth jets - " + range_str
    description = "Filled with the pT response of double-matched reco jets of fCPV < 1 - " + range_str
    xt = "p{T} response of leading jet"
    yt = "Fraction of events"
    plot(group[7], True, title, ROOT.kBlack, False, False, xt, yt, description)


# PLOT PROBABILITY THAT ANY MATCHED JET IS DOUBLE MATCHED ACCORDING TO pT

# Any fCPV

# according to reco pT
efficiency = ROOT.TEfficiency(hist_double_any, hist_all_any)
description = "Filled with fraction of: no. of double matched jets / total no. of matched jets of any fCPV; for each reco pT bin"
title = "Probability of a matched jet being double matched for any fCPV; reco pT; Probability"
plot(efficiency, False, title, ROOT.kRed, False, True, "", "", description)

# according to truth pT
efficiency = ROOT.TEfficiency(hist_double_any_t, hist_all_any_t)
description = "Filled with fraction of: no. of double matched jets / total no. of matched jets of any fCPV; for each total truth-matched pT bin"
title = "Probability of a matched jet being double matched for any fCPV; total truth-matched pT; Probability"
plot(efficiency, False, title, ROOT.kRed, False, True, "", "", description)


# fCPV < 1

# according to reco pT
efficiency = ROOT.TEfficiency(hist_double, hist_all)
description = "Filled with fraction of: no. of double matched jets / total no. of matched jets of fCPV < 1; for each reco pT bin"
title = "Probability of a matched jet being double matched for fCPV < 1; reco pT; Probability"
plot(efficiency, False, title, ROOT.kBlue, False, True, "", "", description)

# according to truth pT
efficiency = ROOT.TEfficiency(hist_double_t, hist_all_t)
description = "Filled with fraction of: no. of double matched jets / total no. of matched jets of fCPV < 1; for each total truth-matched pT bin"
title = "Probability of a matched jet being double matched for fCPV < 1; total truth-matched pT; Probability"
plot(efficiency, False, title, ROOT.kBlue, False, True, "", "", description)

# fCPV = 1

# according to reco pT
efficiency = ROOT.TEfficiency(hist_double_1, hist_all_1)
description = "Filled with fraction of: no. of double matched jets / total no. of matched jets of fCPV = 1; for each reco pT bin"
title = "Probability of a matched jet being double matched for fCPV = 1; reco pT; Probability"
plot(efficiency, False, title, ROOT.kGreen, False, True, "", "", description)

# according to truth pT
efficiency = ROOT.TEfficiency(hist_double_1_t, hist_all_1_t)
description = "Filled with fraction of: no. of double matched jets / total no. of matched jets of fCPV = 1; for each total truth-matched pT bin"
title = "Probability of a matched jet being double matched for fCPV = 1; total truth-matched pT; Probability"
plot(efficiency, False, title, ROOT.kGreen, False, True, "", "", description)


# SAVE PLOTS TO PDF
canvas.Print(outFileName+"]")

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
canvas = ROOT.TCanvas("canvas","canvas",1080,650)
canvas.cd()
canvas.Print(outFileName+"[")
text_size = 0.05

# PLOT PROBABILITY THAT THERE IS NO JET WITHIN dR < 0.4 OF A JET OF fCPV < 1
efficiency = ROOT.TEfficiency(hist_num_deltaR, hist_all_fCPV)
efficiency.SetTitle("Frequency of there being no other jets within dR < 0.4 for a reco jet of fCPV < 1; reco pT; Probability")
canvas.SetLogx(False)
canvas.SetLogy(False)
efficiency.SetLineColor(ROOT.kBlack)
efficiency.Draw("AP")
canvas.Print(outFileName)

# PLOT LEADING JET VS SUBLEADING JET FOR DOUBLE MATCHED JETS

# Any fCPV
canvas.SetLogx(False)
canvas.SetLogy(False)
hist_pt_comparison_fCPV_any.GetXaxis().SetTitle("p{T} of Leading Jet")
hist_pt_comparison_fCPV_any.GetYaxis().SetTitle("p{T} of Subleading Jet")
hist_pt_comparison_fCPV_any.SetTitle("Correlation between Leading Jet p_{T} and Subleading Jet p_{T} for any fCPV")
hist_pt_comparison_fCPV_any.SetStats(0)
hist_pt_comparison_fCPV_any.Draw("colz")
canvas.Print(outFileName)

# fCPV = 1
canvas.SetLogx(False)
canvas.SetLogy(False)
hist_pt_comparison_fCPV_1.GetXaxis().SetTitle("p{T} of Leading Jet")
hist_pt_comparison_fCPV_1.GetYaxis().SetTitle("p{T} of Subleading Jet")
hist_pt_comparison_fCPV_1.SetTitle("Correlation between Leading Jet p_{T} and Subleading Jet p_{T} for fCPV = 1")
hist_pt_comparison_fCPV_1.SetStats(0)
hist_pt_comparison_fCPV_1.Draw("colz")
canvas.Print(outFileName)

# fCPV < 1
canvas.SetLogx(False)
canvas.SetLogy(False)
hist_pt_comparison_fCPV.GetXaxis().SetTitle("p{T} of Leading Jet")
hist_pt_comparison_fCPV.GetYaxis().SetTitle("p{T} of Subleading Jet")
hist_pt_comparison_fCPV.SetTitle("Correlation between Leading Jet p_{T} and Subleading Jet p_{T} for fCPV < 1")
hist_pt_comparison_fCPV.SetStats(0)
hist_pt_comparison_fCPV.Draw("colz")
canvas.Print(outFileName)


# PLOT PROBABILITY THAT ANY MATCHED JET IS DOUBLE MATCHED ACCORDING TO pT
efficiency_any = ROOT.TEfficiency(hist_double_any, hist_all_any)
efficiency = ROOT.TEfficiency(hist_double, hist_all)
efficiency_1 = ROOT.TEfficiency(hist_double_1, hist_all_1)

efficiency_any_t = ROOT.TEfficiency(hist_double_any_t, hist_all_any_t)
efficiency_t = ROOT.TEfficiency(hist_double_t, hist_all_t)
efficiency_1_t = ROOT.TEfficiency(hist_double_1_t, hist_all_1_t)

efficiency_any.SetTitle("Probability of a matched jet being double matched for any fCPV; reco pT; Probability")
efficiency.SetTitle("Probability of a matched jet being double matched for fCPV < 1; reco pT; Probability")
efficiency_1.SetTitle("Probability of a matched jet being double matched for fCPV = 1; reco pT; Probability")

efficiency_any_t.SetTitle("Probability of a matched jet being double matched for any fCPV; total truth-matched pT; Probability")
efficiency_t.SetTitle("Probability of a matched jet being double matched for fCPV < 1; total truth-matched pT; Probability")
efficiency_1_t.SetTitle("Probability of a matched jet being double matched for fCPV = 1; total truth-matched pT; Probability")

# Any fCPV

# according to reco pT
canvas.SetLogx(False)
canvas.SetLogy(False)
efficiency_any.SetLineColor(ROOT.kRed)
efficiency_any.Draw("AP")
canvas.Print(outFileName)

# according to truth pT
canvas.SetLogx(False)
canvas.SetLogy(False)
efficiency_any_t.SetLineColor(ROOT.kRed)
efficiency_any_t.Draw("AP")
canvas.Print(outFileName)

# fCPV < 1

# according to reco pT
canvas.SetLogx(False)
canvas.SetLogy(False)
efficiency.SetLineColor(ROOT.kBlue)
efficiency.Draw("AP")
canvas.Print(outFileName)

# according to truth pT
canvas.SetLogx(False) 
canvas.SetLogy(False)
efficiency_t.SetLineColor(ROOT.kBlue)
efficiency_t.Draw("AP")
canvas.Print(outFileName) 

# fCPV = 1

# according to reco pT
canvas.SetLogx(False)
canvas.SetLogy(False)
efficiency_1.SetLineColor(ROOT.kGreen)
efficiency_1.Draw("AP")
canvas.Print(outFileName)

# according to truth pT
canvas.SetLogx(False)
canvas.SetLogy(False)
efficiency_1_t.SetLineColor(ROOT.kGreen)
efficiency_1_t.Draw("AP")
canvas.Print(outFileName) 



# PLOT pT RESPONSE PLOTS
for j in range(len(ranges)):
    group = hists_pt_responses[j]
    range_str = "Range: " + str(ranges[j][0]) + "-" + str(ranges[j][1]) + " MeV"

    either_or = group[0].Clone("either_or")
    either_or.Add(group[3])

    # ANY fCPV

    # only truth-matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    group[0].Scale(1./group[0].Integral())
    group[0].SetLineColor(ROOT.kBlack)
    group[0].SetLineWidth(2)
    group[0].SetTitle("pT response: p{T} of reco jets of any fCPV / p{T} of matching truth jet - " + range_str)
    group[0].GetXaxis().SetTitle("p{T} response of leading jet")
    group[0].GetYaxis().SetTitle("Fraction of events")
    group[0].SetStats(0)
    fit_gaus(group[0])
    canvas.Print(outFileName)

    # only pileup-matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    group[3].Scale(1./group[3].Integral())
    group[3].SetLineColor(ROOT.kBlack)
    group[3].SetLineWidth(2)
    group[3].SetTitle("p{T} response: p{T} of reco jets of any fCPV / p{T} of matching inTimeTruthJet - " + range_str)
    group[3].GetXaxis().SetTitle("p{T} response of leading jet")
    group[3].GetYaxis().SetTitle("Fraction of events")
    group[3].SetStats(0)
    fit_gaus(group[3])
    canvas.Print(outFileName)

    # matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    either_or.Scale(1./either_or.Integral())
    either_or.SetLineColor(ROOT.kBlack)
    either_or.SetLineWidth(2)
    either_or.SetTitle("p{T} response: p{T} of reco jets of any fCPV / p{T} of matching truth jet OR inTimeTruthJet - " + range_str)
    either_or.GetXaxis().SetTitle("p{T} response of leading jet")
    either_or.GetYaxis().SetTitle("Fraction of events")
    either_or.SetStats(0)
    fit_gaus(either_or)
    canvas.Print(outFileName)

    # double matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    group[6].Scale(1./group[6].Integral())
    group[6].SetLineColor(ROOT.kBlack)
    group[6].SetLineWidth(2)
    group[6].SetTitle("p{T} response: p{T} of reco jets of any fCPV / p{T} sum of BOTH matching truth jets - " + range_str)
    group[6].GetXaxis().SetTitle("p{T} response of leading jet")
    group[6].GetYaxis().SetTitle("Fraction of events")
    group[6].SetStats(0)
    fit_gaus(group[6])
    canvas.Print(outFileName)

    # fCPV = 1
    either_or = group[2].Clone("either_or")
    either_or.Add(group[5])

    # only truth-matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    group[2].Scale(1./group[2].Integral())
    group[2].SetLineColor(ROOT.kBlack)
    group[2].SetLineWidth(2)
    group[2].SetTitle("p{T} response: p{T} of reco jets of fCPV = 1 / p{T} of matching truth jet - " + range_str)
    group[2].GetXaxis().SetTitle("p{T} response of leading jet")
    group[2].GetYaxis().SetTitle("Fraction of events")
    group[2].SetStats(0)
    fit_gaus(group[2])
    canvas.Print(outFileName)

    # only pileup-matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    group[5].Scale(1./group[5].Integral())
    group[5].SetLineColor(ROOT.kBlack)
    group[5].SetLineWidth(2)
    group[5].SetTitle("p{T} response: p{T} of reco jets of fCPV = 1 / p{T} of matching inTimeTruthJet - " + range_str)
    group[5].GetXaxis().SetTitle("p{T} response of leading jet")
    group[5].GetYaxis().SetTitle("Fraction of events")
    group[5].SetStats(0)
    fit_gaus(group[5])
    canvas.Print(outFileName)

    # matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    either_or.Scale(1./either_or.Integral())
    either_or.SetLineColor(ROOT.kBlack)
    either_or.SetLineWidth(2)
    either_or.SetTitle("p{T} response: p{T} of reco jets of fCPV = 1 / p{T} of matching truth jet OR inTimeTruthJet - " + range_str)
    either_or.GetXaxis().SetTitle("p{T} response of leading jet")
    either_or.GetYaxis().SetTitle("Fraction of events")
    either_or.SetStats(0)
    fit_gaus(either_or)
    canvas.Print(outFileName)

    # double matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    group[8].Scale(1./group[8].Integral())
    group[8].SetLineColor(ROOT.kBlack)
    group[8].SetLineWidth(2)
    group[8].SetTitle("p{T} response: p{T} of reco jets of fCPV = 1 / p{T} sum of BOTH matching truth jets - " + range_str)
    group[8].GetXaxis().SetTitle("p{T} response of leading jet")
    group[8].GetYaxis().SetTitle("Fraction of events")
    group[8].SetStats(0)
    fit_gaus(group[8])
    canvas.Print(outFileName)

    # fCPV < 1
    either_or = group[1].Clone("either_or")
    either_or.Add(group[4])

    # only truth-matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    group[1].Scale(1./group[1].Integral())
    group[1].SetLineColor(ROOT.kBlack)
    group[1].SetLineWidth(2)
    group[1].SetTitle("p{T} response: p{T} of leading reco jet of fCPV < 1 / p{T} of matching truth jet - " + range_str)
    group[1].GetXaxis().SetTitle("p{T} response of leading jet")
    group[1].GetYaxis().SetTitle("Fraction of events")
    group[1].SetStats(0)
    fit_gaus(group[1])
    canvas.Print(outFileName)

    # only pileup-matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    group[4].Scale(1./group[4].Integral())
    group[4].SetLineColor(ROOT.kBlack)
    group[4].SetLineWidth(2)
    group[4].SetTitle("p{T} response: p{T} of leading reco jet of fCPV < 1 / p{T} of matching inTimeTruthJet - " + range_str)
    group[4].GetXaxis().SetTitle("p{T} response of leading jet")
    group[4].GetYaxis().SetTitle("Fraction of events")
    group[4].SetStats(0)
    fit_gaus(group[4])
    canvas.Print(outFileName)

    # matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    either_or.Scale(1./either_or.Integral())
    either_or.SetLineColor(ROOT.kBlack)
    either_or.SetLineWidth(2)
    either_or.SetTitle("p{T} response: p{T} of leading reco jet of fCPV < 1 / p{T} of matching truth jet OR inTimeTruthJet - " + range_str)
    either_or.GetXaxis().SetTitle("p{T} response of leading jet")
    either_or.GetYaxis().SetTitle("Fraction of events")
    either_or.SetStats(0)
    fit_gaus(either_or)
    canvas.Print(outFileName)

    # double matched
    canvas.SetLogx(False)
    canvas.SetLogy(False)
    group[7].Scale(1./group[7].Integral())
    group[7].SetLineColor(ROOT.kBlack)
    group[7].SetLineWidth(2)
    group[7].SetTitle("p{T} response: p{T} of leading reco jet of fCPV < 1 / p{T} sum of BOTH matching truth jets - " + range_str)
    group[7].GetXaxis().SetTitle("p{T} response of leading jet")
    group[7].GetYaxis().SetTitle("Fraction of events")
    group[7].SetStats(0)
    fit_gaus(group[7])
    canvas.Print(outFileName)



# PLOT FCPV COMPARISON HISTOGRAM
canvas.SetLogx(False)
canvas.SetLogy(False)
hist_cFPV_unf.SetLineColor(ROOT.kRed)
hist_cFPV_unf.SetLineWidth(2)
hist_cFPV_unf.GetXaxis().SetTitle("fCPV of reco jets")
hist_cFPV_unf.GetYaxis().SetTitle("Scaled number of events")
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
hist_fCPV_truth_pileup_or.GetYaxis().SetTitle("Scaled number of events")
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
canvas.Print(outFileName)

canvas.Print(outFileName+"]")

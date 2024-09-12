import ROOT
import fastjet
import sys
import math
import numpy as np

# Load input files and TTrees
inFile_names = [

	#"data/group.phys-exotics.41005746.EXT0._000001.ByVertexJetNtuple.root",
	"data/group.phys-exotics.41005731.EXT0._000001.ByVertexJetNtuple.root",
	"data/group.phys-exotics.41005731.EXT0._000002.ByVertexJetNtuple.root",
	"data/group.phys-exotics.41005746.EXT0._000002.ByVertexJetNtuple.root"

	]


inTree = ROOT.TChain("Vertex")
for inFile_name in inFile_names:
    inTree.Add(inFile_name)

inTree.SetBranchStatus("*", 1) # Enable reading of branches

# Initialise histograms
hist_fCPV_unf = ROOT.TH1F("fCPV_unfiltered","fCPV of all reconstructed jets in all vertices",50,0,1.1)
hist_fCPV_truth = ROOT.TH1F("fCPV_truth","fCPV of all truth-matched reconstructed jets in all vertices",50,0,1.1)
hist_fCPV_pileup = ROOT.TH1F("fCPV_pileup","fCPV of all truth-pileup-matched reconstructed jets in all vertices",50,0,1.1)
hist_fCPV_truth_pileup = ROOT.TH1F("fCPV_truth_pileup","fCPV of all truth- and pileup-matched reconstructed jets in all vertices",50,0,1.1)
hist_fCPV_truth_pileup_or = ROOT.TH1F("fCPV_truth_pileup_or","fCPV of all truth- or pileup-matched reconstructed jets in all vertices",50,0,1.1)

hist_pt_comparison_any = ROOT.TH2I("pt_comparison_any","Correlation between jet 1 p_{T} and jet 2 p_{T}",75,10.e3,40.e3,75,10.e3,40.e3)
hist_pt_comparison = ROOT.TH2I("pt_comparison","Correlation between jet 1 p_{T} and jet 2 p_{T}",75,10.e3,40.e3,75,10.e3,40.e3)
hist_pt_comparison_1 = ROOT.TH2I("pt_comparison_1","Correlation between jet 1 p_{T} and jet 2 p_{T}",75,10.e3,40.e3,75,10.e3,40.e3)

hist_double_any = ROOT.TH1F("double_any","Number of jets in pT ranges that are double matched for any fCPV",50,10.e3,40.e3)
hist_all_any = ROOT.TH1F("all_any","Number of jets in pT ranges for any fCPV",50,10.e3,40.e3)

hist_double = ROOT.TH1F("double","Number of jets in pT ranges that are double matched for fCPV < 1",50,10.e3,40.e3)
hist_all = ROOT.TH1F("all","Number of jets in pT ranges for fCPV < 1",50,10.e3,40.e3)

hist_double_1 = ROOT.TH1F("double_1","Number of jets in pT ranges that are double matched for fCPV = 1",50,10.e3,40.e3)
hist_all_1 = ROOT.TH1F("all_1","Number of jets in pT ranges for fCPV = 1",50,10.e3,40.e3)

hist_double_any_t = ROOT.TH1F("double_any_t","Number of jets in pT ranges that are double matched for any fCPV",50,10.e3,40.e3)
hist_all_any_t = ROOT.TH1F("all_any_t","Number of jets in pT ranges for any fCPV",50,10.e3,40.e3)

hist_double_t = ROOT.TH1F("double_t","Number of jets in pT ranges that are double matched for fCPV < 1",50,10.e3,40.e3)
hist_all_t = ROOT.TH1F("all_t","Number of jets in pT ranges for fCPV < 1",50,10.e3,40.e3)

hist_double_1_t = ROOT.TH1F("double_1_t","Number of jets in pT ranges that are double matched for fCPV = 1",50,10.e3,40.e3)
hist_all_1_t = ROOT.TH1F("all_1_t","Number of jets in pT ranges for fCPV = 1",50,10.e3,40.e3)

ranges = [[20.e3, 25.e3], [25.e3, 30.e3], [30.e3, 40.e3]]
hists_pt_responses = []
for i in ranges:
	base = str(i[0]) + "-" + str(i[1])
	range_str = "Range = " + base + " MeV"
	for_range = []
	for_range.append(ROOT.TH1F("pt_response_any_"+base,"pt response of leading reconstructed jets from vertices; " + range_str,50,0,2))
	for_range.append(ROOT.TH1F("pt_response_"+base,"pt response of leading reconstructed jets from vertices; " + range_str,50,0,2))
	for_range.append(ROOT.TH1F("pt_response_1_"+base,"pt response of leading reconstructed jets (with fCPV) from vertices; " + range_str,50,0,2))

	for_range.append(ROOT.TH1F("pt_response_any_p_"+base,"pt response of leading reconstructed jets from vertices; " + range_str,50,0,2))
	for_range.append(ROOT.TH1F("pt_response_p_"+base,"pt response of leading reconstructed jets from vertices; " + range_str,50,0,2))
	for_range.append(ROOT.TH1F("pt_response_1_p_"+base,"pt response of leading reconstructed jets (with fCPV) from vertices; " + range_str,50,0,2))

	for_range.append(ROOT.TH1F("pt_response_any_both_"+base,"pt response of leading reconstructed jets from vertices; " + range_str,50,0,2))
	for_range.append(ROOT.TH1F("pt_response_both_"+base,"pt response of leading reconstructed jets from vertices; " + range_str,50,0,2))
	for_range.append(ROOT.TH1F("pt_response_1_both_"+base,"pt response of leading reconstructed jets (with fCPV) from vertices; " + range_str,50,0,2))

	hists_pt_responses.append(for_range)


# Calculates dR as a function of dEta and dPhi
def calc_delta_r(phi1, phi2, eta1, eta2):
	delta_phi = phi2-phi1
	delta_eta = eta2-eta1
	return math.sqrt(delta_phi**2 + delta_eta**2)


jets_temp = []
prev = -1

num_entries = inTree.GetEntries()
temp = 1000000

# Process entries
for i in range(num_entries):
	if i%2000 == 0:
		print(f"{i} / {num_entries}")

	inTree.GetEntry(i)

	# Process the truth and reco jets in each vertex
	for j in range(inTree.jet_truthPt.size()):
		jet_pt_frac = inTree.jet_ptFrac_dR04.at(j) # = fCPV

		truthJet_pt = inTree.jet_truthPt.at(j)
		truthPileupJet_pt = inTree.jet_truthPileupPt.at(j)

		total_pt = 0
		if truthJet_pt != -999.0:
			total_pt += truthJet_pt
		if truthPileupJet_pt != -999.0:
			total_pt += truthPileupJet_pt

		is_truth_matched = np.bool_(inTree.jet_isTruthMatched.at(j))
		is_pileup_matched = np.bool_(inTree.jet_isTruthPileupMatched.at(j))

		is_truth_matched = is_truth_matched and (inTree.jet_truthPt.at(j) > 10.e3)
		is_pileup_matched = is_pileup_matched and (inTree.jet_truthPileupPt.at(j) > 10.e3)
		double_matched = is_truth_matched and is_pileup_matched
		matched = is_truth_matched or is_pileup_matched
		
		# Fill pT response histograms according to: a) what truth jet(s) the reco jet is matched to; b) the pT of the matched truth jet(s)
		if matched:
			# Calculate pT responses
			ratio = inTree.jet_pt.at(j) / truthJet_pt # only truth matched
			ratio_p = inTree.jet_pt.at(j) / truthPileupJet_pt # only pileup matched
			ratio_both = inTree.jet_pt.at(j) / total_pt # double matched

			# Fill the appropriate pT range histogram depending on the pT of matched jet(s)
			for k in range(len(ranges)): # 20-25 GeV; 25-30 GeV; 30-40 GeV
				# Check if the jet fits the pT range
				pt_filter = ranges[k][0] < total_pt < ranges[k][1]
				if pt_filter:
					# Find appropriate histogram to fill, for any fCPV, and then within fCPV = 1 OR fCPV < 1
					if is_truth_matched and not is_pileup_matched:
						hists_pt_responses[k][0].Fill(ratio)
					if is_pileup_matched and not is_truth_matched:
						hists_pt_responses[k][3].Fill(ratio_p)
					if is_pileup_matched and is_truth_matched:
						hists_pt_responses[k][6].Fill(ratio_both)

					if jet_pt_frac == 1: # fCPV = 1
						if is_truth_matched and not is_pileup_matched:
							hists_pt_responses[k][2].Fill(ratio)
						if is_pileup_matched and not is_truth_matched:
							hists_pt_responses[k][5].Fill(ratio_p)
						if is_pileup_matched and is_truth_matched:
							hists_pt_responses[k][8].Fill(ratio_both)

					elif jet_pt_frac < 1: # fCPV < 1
						if is_truth_matched and not is_pileup_matched:
							hists_pt_responses[k][1].Fill(ratio)
						if is_pileup_matched and not is_truth_matched:
							hists_pt_responses[k][4].Fill(ratio_p)
						if is_pileup_matched and is_truth_matched:
							hists_pt_responses[k][7].Fill(ratio_both)

		# Fill the 1D fCPV histograms for: truth-matched, pileup-matched, and unfiltered
		hist_fCPV_unf.Fill(jet_pt_frac)
		if matched:
			hist_fCPV_truth_pileup_or.Fill(jet_pt_frac)

		if double_matched:
			hist_fCPV_truth_pileup.Fill(jet_pt_frac)

		elif is_truth_matched and not is_pileup_matched:
			hist_fCPV_truth.Fill(jet_pt_frac)

		elif is_pileup_matched and not is_truth_matched:
			hist_fCPV_pileup.Fill(jet_pt_frac)

		

		# Fill 1D pT histograms for matched truth jet(s) and reco jet
		if matched:
			hist_all_any.Fill(inTree.jet_pt.at(j)) # reco jet
			hist_all_any_t.Fill(total_pt) # matched truth jet(s)

		# Fill 1D pT histograms for matched truth jets and reco jet which is DOUBLE matched
		if double_matched:
			hist_double_any.Fill(inTree.jet_pt.at(j))
			hist_double_any_t.Fill(total_pt)

			# Also, fill 2D pT histogram of leading and subleading jet pTs of both truth jets that work to form the reco jet (also for fCPV = 1, < 1, and inclusive)
			leading = max(truthJet_pt, truthPileupJet_pt)
			sleading = min(truthJet_pt, truthPileupJet_pt)

			hist_pt_comparison_any.Fill(leading, sleading)
			if (jet_pt_frac == 1):
				hist_pt_comparison_1.Fill(leading, sleading)
			elif (jet_pt_frac < 1):
				hist_pt_comparison.Fill(leading, sleading)

		# Fill 1D pT histograms for matched truth jet(s) and reco jet as well as for the jets that are double matched, for fCPV = 1 and fCPV < 1
		if jet_pt_frac < 1:
			if matched:
				hist_all.Fill(inTree.jet_pt.at(j))
				hist_all_t.Fill(total_pt)
			if (double_matched):
				hist_double.Fill(inTree.jet_pt.at(j))
				hist_double_t.Fill(total_pt)

		elif jet_pt_frac == 1:
			if matched:
				hist_all_1.Fill(inTree.jet_pt.at(j))
				hist_all_1_t.Fill(total_pt)
			if (double_matched):
				hist_double_1.Fill(inTree.jet_pt.at(j))
				hist_double_1_t.Fill(total_pt)


# Write all histograms to root file
outFile = ROOT.TFile.Open("ntuple_analysis.root","RECREATE")
outFile.cd()

hist_fCPV_unf.Write()
hist_fCPV_truth.Write()
hist_fCPV_pileup.Write()
hist_fCPV_truth_pileup.Write()
hist_fCPV_truth_pileup_or.Write()

hist_pt_comparison.Write()
hist_pt_comparison_1.Write()
hist_pt_comparison_any.Write()

for j in hists_pt_responses:
	for k in j:
		k.Write()

hist_double_any.Write()
hist_all_any.Write()

hist_double.Write()
hist_all.Write()

hist_double_1.Write()
hist_all_1.Write()

hist_double_any_t.Write()
hist_all_any_t.Write()

hist_double_t.Write()
hist_all_t.Write()

hist_double_1_t.Write()
hist_all_1_t.Write()

outFile.Close()

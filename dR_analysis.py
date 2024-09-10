import ROOT
import fastjet
import sys
import math
import numpy as np
from collections import defaultdict

# Load input files and TTree
inFile = ROOT.TFile.Open("user.cmorenom.41078044.EXT0._000001.TruthJetVertexingNtuple.root","READ")

inTree = inFile.Get("Event")
if not inTree:
    printf("Failed to retrieve the input tree from input file:",inFileName)
    sys.exit(1)

# Initialise histograms
hist_deltaR_04 = ROOT.TH1F("deltaR04","dR of overlapping truth jets",200,0,0.5)
hist_pt_deltaR = ROOT.TH2I("deltaR_pt","Correlation between jet 1 p_{T} and jet 2 p_{T}",75,0,0.4,75,1,10)
hist_all_pt = ROOT.TH1F("all_pt","pT of all truth jets",200,10.e3,40.e3)
hist_dR_pt = ROOT.TH1F("dR_pt","pT of all overlapping truth jets",200,10.e3,40.e3)

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

hist_num_deltaR = ROOT.TH1F("num_deltaR","Number of jets in pT ranges for fCPV < 1 and dR < 0.4 with any neighbouring jet",50,10.e3,40.e3)
hist_all_fCPV = ROOT.TH1F("all_fCPV","Number of jets in pT ranges for fCPV < 1",50,10.e3,40.e3)

# Enable reading of branches
inTree.SetBranchStatus("*", 1)


# Calculates dR as a function of dEta and dPhi
def calc_delta_r(phi1, phi2, eta1, eta2):
	delta_phi = phi2-phi1
	delta_eta = eta2-eta1
	return math.sqrt(delta_phi**2 + delta_eta**2)


# Process entries
temp = 100
num_entries = inTree.GetEntries()
for i in range(num_entries): # 'temp' for testing; 'num_entries' for real data
	if i%200 == 0:
		print(f"{i} / {num_entries}")
	
	inTree.GetEntry(i)
	truth_jets = []

	# Loop through trutJets and truthInTimeJets and add save their basic properties in a list
	for j in range(inTree.truthJet_pt.size()):
		jet = [inTree.truthJet_pt.at(j), inTree.truthJet_eta.at(j), inTree.truthJet_phi.at(j)]
		truth_jets.append(jet)

	for j in range(inTree.truthInTimeJet_pt.size()):
		jet = [inTree.truthInTimeJet_pt.at(j), inTree.truthInTimeJet_eta.at(j), inTree.truthInTimeJet_phi.at(j)]
		truth_jets.append(jet)

	# Calculate dR between all truth jets
	dR_matrix = [[999 for _ in range(len(truth_jets))] for _ in range(len(truth_jets))]
	for j in range(len(truth_jets)):
		hist_all_pt.Fill(truth_jets[j][0])
		eta1, phi1 = truth_jets[j][1], truth_jets[j][2]
		for k in range(j+1, len(truth_jets)):
			eta2, phi2 = truth_jets[k][1], truth_jets[k][2]
			delta_r = calc_delta_r(phi1, phi2, eta1, eta2)
			dR_matrix[j][k] = delta_r
			if delta_r < 0.4:
				ratio = max(truth_jets[j][0], truth_jets[k][0]) / min(truth_jets[j][0], truth_jets[k][0])
				hist_deltaR_04.Fill(delta_r)

		# a jet is within dR < 0.4 of any other jet, it overlaps, therefore fill appropriate histogram with pT
		if np.any(np.array(dR_matrix[j]) < 0.4):
			hist_dR_pt.Fill(truth_jets[j][0])

	# Find jets with fCPV < 1, fill list with all their properties
	jets = []
	for j in range(inTree.jet_pt.size()):
		if inTree.jet_ptFrac_dR04.at(j) < 1:
			jet = [
				inTree.jet_pt.at(j),
				inTree.jet_eta.at(j),
				inTree.jet_phi.at(j),
				inTree.jet_originVertexIndex.at(j),
				]
			jets.append(jet)
			hist_all_fCPV.Fill(jet[0])

	# Calculate dR between all these reco (fCPV < 1) jets
	dR_matrix_reco = [[999 for _ in range(len(jets))] for _ in range(len(jets))]
	jets_properties = np.array(jets)
	phis = jets_properties[:, 2]
	etas = jets_properties[:, 1]

	zipped_list = [[float(e), float(p)] for e, p in zip(etas, phis)]

	for j in range(len(jets)):
		for k in range(0, len(jets)): #0 or j+1
			if j == k or jets[j][3] == jets[k][3]:
				continue

			eta1, phi1 = zipped_list[j]
			eta2, phi2 = zipped_list[k]
			delta_r = calc_delta_r(phi1, phi2, eta1, eta2)
			dR_matrix_reco[j][k] = delta_r

	# if there is no other jet within dR < 0.4, fill appropriate histogram with pT
	for j in range(len(dR_matrix_reco)):
		pt = jets[j][0]
		count = sum([1 for x in dR_matrix_reco[j] if x < 0.4])
		if count == 0:
			hist_num_deltaR.Fill(pt)




# Write histograms to file
outFile = ROOT.TFile.Open("truthDR_analysis.root","RECREATE")
outFile.cd()

hist_deltaR_04.Write()
hist_pt_deltaR.Write()
hist_all_pt.Write()
hist_dR_pt.Write()

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

hist_num_deltaR.Write()
hist_all_fCPV.Write()

outFile.Close()

import ROOT
import fastjet
import sys
import math
import numpy as np
from collections import defaultdict

# Load input files and TTree
inFile_names = [

	"data/user.vcepaiti.41131452.EXT0._000001.TruthJetVertexingNtuple.root",
	"data/user.vcepaiti.41131452.EXT0._000002.TruthJetVertexingNtuple.root",
	"data/user.vcepaiti.41131452.EXT0._000003.TruthJetVertexingNtuple.root",
	#"data/user.vcepaiti.41131452.EXT0._000004.TruthJetVertexingNtuple.root",
	#"data/user.vcepaiti.41131452.EXT0._000005.TruthJetVertexingNtuple.root",
	#"data/user.vcepaiti.41131452.EXT0._000006.TruthJetVertexingNtuple.root",
	#"data/user.vcepaiti.41131452.EXT0._000007.TruthJetVertexingNtuple.root",

	]


inTree = ROOT.TChain("Event")
for inFile_name in inFile_names:
    inTree.Add(inFile_name)


# Initialise histograms
hist_deltaR_04 = ROOT.TH1F("deltaR04","dR of overlapping truth jets",200,0,0.5)
hist_pt_deltaR = ROOT.TH2I("deltaR_pt","Correlation between jet 1 p_{T} and jet 2 p_{T}",75,0,0.4,75,1,10)

hist_all_pt = ROOT.TH1F("all_pt","pT of all truth jets",200,10.e3,40.e3)
hist_dR_pt = ROOT.TH1F("dR_pt","pT of all overlapping truth jets",200,10.e3,40.e3)

hist_double = ROOT.TH1F("double","Number of jets in pT ranges that are double matched for fCPV < 1",50,10.e3,40.e3)
hist_double_1 = ROOT.TH1F("double_1","Number of jets in pT ranges that are double matched for fCPV = 1",50,10.e3,40.e3)
hist_double_p = ROOT.TH1F("double_p","Number of jets in pT ranges that are double pileup-matched for fCPV < 1",50,10.e3,40.e3)
hist_double_1_p = ROOT.TH1F("double_1_p","Number of jets in pT ranges that are double pileup-matched for fCPV < 1",50,10.e3,40.e3)

hist_gpv = ROOT.TH1F("gpv","Number of jets in pT ranges that are double matched for any fCPV",50,10.e3,40.e3)
hist_pileup = ROOT.TH1F("pileup","Number of jets in pT ranges for any fCPV",50,10.e3,40.e3)

hist_gpv_1 = ROOT.TH1F("gpv_1","Number of jets in pT ranges that are double matched for any fCPV",50,10.e3,40.e3)
hist_pileup_1 = ROOT.TH1F("pileup_1","Number of jets in pT ranges for any fCPV",50,10.e3,40.e3)

hist_num_deltaR = ROOT.TH1F("num_deltaR","Number of jets in pT ranges for fCPV < 1 and dR < 0.4 with any neighbouring jet",50,10.e3,40.e3)
hist_all_fCPV = ROOT.TH1F("all_fCPV","Number of jets in pT ranges for fCPV < 1",50,10.e3,40.e3)

hist_num_deltaR_1 = ROOT.TH1F("num_deltaR_1","Number of jets in pT ranges for fCPV = 1 and dR < 0.4 with any neighbouring jet",50,10.e3,40.e3)
hist_all_fCPV_1 = ROOT.TH1F("all_fCPV_1","Number of jets in pT ranges for fCPV = 1",50,10.e3,40.e3)

hist_all_filt = ROOT.TH1F("all_fCPV_filt","Number of jets in pT ranges for fCPV < 1",50,10.e3,40.e3)
hist_all_filt_1 = ROOT.TH1F("all_fCPV_filt_1","Number of jets in pT ranges for fCPV = 1",50,10.e3,40.e3)

# Enable reading of branches
inTree.SetBranchStatus("*", 1)


# Calculates dR as a function of dEta and dPhi
def calc_delta_r(pt1, pt2, eta1, eta2, phi1, phi2, m1, m2):
	# Define root four-vectors
	vec1 = ROOT.TLorentzVector()
	vec1.SetPtEtaPhiM(pt1, eta1, phi1, m1)

	vec2 = ROOT.TLorentzVector()
	vec2.SetPtEtaPhiM(pt2, eta2, phi2, m2)

	# Calculate deltaR
	dR = vec1.DeltaR(vec2)
	return dR


# Process entries
temp = 1000
num_entries = inTree.GetEntries()
for i in range(num_entries): # 'temp' for testing; 'num_entries' for real data
	if i%2000 == 0:
		print(f"{i} / {num_entries}")
	
	inTree.GetEntry(i)

	weight = inTree.event_mcEventWeight

	truth_jets = []

	# Loop through trutJets and truthInTimeJets and add save their basic properties in a list
	for j in range(inTree.truthJet_pt.size()):
		jet = [
			inTree.truthJet_pt.at(j),
			inTree.truthJet_eta.at(j),
			inTree.truthJet_phi.at(j)
		]
		truth_jets.append(jet)

	for j in range(inTree.truthInTimeJet_pt.size()):
		jet = [
			inTree.truthInTimeJet_pt.at(j),
			inTree.truthInTimeJet_eta.at(j),
			inTree.truthInTimeJet_phi.at(j)
		]
		truth_jets.append(jet)

	# Calculate dR between all truth jets
	dR_matrix = [[999 for _ in range(len(truth_jets))] for _ in range(len(truth_jets))]
	for j in range(len(truth_jets)):
		hist_all_pt.Fill(truth_jets[j][0], weight)
		pt1, eta1, phi1 = truth_jets[j][0], truth_jets[j][1], truth_jets[j][2]
		for k in range(j+1, len(truth_jets)):
			pt2, eta2, phi2 = truth_jets[k][0], truth_jets[k][1], truth_jets[k][2]
			delta_r = calc_delta_r(pt1, pt2, eta1, eta2, phi1, phi2, 0, 0)
			dR_matrix[j][k] = delta_r
			if delta_r < 0.4:
				ratio = max(truth_jets[j][0], truth_jets[k][0]) / min(truth_jets[j][0], truth_jets[k][0])
				hist_deltaR_04.Fill(delta_r, weight)

		# a jet is within dR < 0.4 of any other jet, it overlaps, therefore fill appropriate histogram with pT
		if np.any(np.array(dR_matrix[j]) < 0.4):
			hist_dR_pt.Fill(truth_jets[j][0], weight)

	# Fill jets[0] with the properties of all jets of fCPV < 1, and same for fCPV = 1 in jets[1]
	jets = [[], []]
	for j in range(inTree.jet_pt.size()):
		jet = [
			inTree.jet_pt.at(j),
			inTree.jet_eta.at(j),
			inTree.jet_phi.at(j),
			inTree.jet_m.at(j),
			inTree.jet_originVertexIndex.at(j),
			inTree.jet_ptFrac_dR04.at(j),
			inTree.jet_Jvt.at(j)
		]
		'''
		control = True
		if inTree.event_truthMatrixElementVertexIndex != inTree.jet_originVertexIndex.at(j):
			jet.append(-1)
			control = False
		'''

		if inTree.jet_ptFrac_dR04.at(j) < 1:
			jets[0].append(jet)
			hist_all_fCPV.Fill(jet[0], weight)
			'''
			if control:
				hist_all_filt.Fill(jet[0], weight)
			'''

		elif inTree.jet_ptFrac_dR04.at(j) == 1:
			jets[1].append(jet)
			hist_all_fCPV_1.Fill(jet[0], weight)
			'''
			if control:
				hist_all_filt_1.Fill(jet[0], weight)
			'''

	# Process all reco jets and calculate dR WRT each other
	for x in range(len(jets)):
		if len(jets[x]) > 0:
			if len(jets[x]) == 1:
				if x == 0:
					hist_num_deltaR.Fill(jets[x][0][0], weight)
				else:
					hist_num_deltaR_1.Fill(jets[x][0][0], weight)

				continue

			# Calculate dR between all these reco (fCPV < 1) jets
			dR_matrix = [[999 for _ in range(len(jets[x]))] for _ in range(len(jets[x]))]

			for j in range(len(jets[x])):
				for k in range(0, len(jets[x])): #0 or j+1
					if j == k or jets[x][j][4] == jets[x][k][4]:
						continue

					pt1, eta1, phi1, m1 = jets[x][j][0:4]
					pt2, eta2, phi2, m2 = jets[x][k][0:4]
					delta_r = calc_delta_r(pt1, pt2, eta1, eta2, phi1, phi2, m1, m2)
					dR_matrix[j][k] = delta_r

			# if there is no other jet within dR < 0.4, fill appropriate histogram with pT
			for j in range(len(dR_matrix)):
				pt = jets[x][j][0]
				count = sum([1 for x in dR_matrix[j] if x < 0.4])
				if count == 0:
					if x == 0: # fCPV < 1
						hist_num_deltaR.Fill(pt, weight)
					else: # fCPV = 1
						hist_num_deltaR_1.Fill(pt, weight)

	# Process all reco jets and calculate dR WRT truth jets, to see if the reco jets match (dR < 0.3)
	for x in range(len(jets)):
		for j in range(len(jets[x])):
			pt1, eta1, phi1, m1 = jets[x][j][0:4]
			indices = []
			for k in range(len(truth_jets)):
				pt2, eta2, phi2 = truth_jets[k]
				delta_r = calc_delta_r(pt1, pt2, eta1, eta2, phi1, phi2, m1, 0)
				if delta_r < 0.3:
					#matches += 1
					indices.append(k)

			gpv_matched = np.any(np.array(indices) < inTree.truthJet_pt.size())
			indices = np.array(indices)
			pileup_matched = indices >= inTree.truthJet_pt.size()
			double_pileup = np.count_nonzero(pileup_matched) >= 2

			# Depending on to what extent the reco jet matches, fill different histograms
			if gpv_matched:
				if x == 0: # fCPV < 1
					hist_gpv.Fill(jets[0][j][0], weight)

				else: # fCPV = 1
					hist_gpv_1.Fill(jets[1][j][0], weight)

			if np.any(pileup_matched):
				if x == 0:
					hist_pileup.Fill(jets[0][j][0], weight)

				else:
					hist_pileup_1.Fill(jets[1][j][0], weight)

			if double_pileup:
				if x == 0:
					hist_double_p.Fill(jets[0][j][0], weight)

				else:
					hist_double_1_p.Fill(jets[1][j][0], weight)

			if gpv_matched and np.any(pileup_matched):
				if x == 0:
					hist_double.Fill(jets[0][j][0], weight)
					
				else:
					hist_double_1.Fill(jets[1][j][0], weight)


# Write histograms to file
outFile = ROOT.TFile.Open("dR_analysis.root","RECREATE")
outFile.cd()

hist_deltaR_04.Write()
hist_pt_deltaR.Write()
hist_all_pt.Write()
hist_dR_pt.Write()

hist_double.Write()
hist_double_1.Write()

hist_num_deltaR.Write()
hist_all_fCPV.Write()

hist_num_deltaR_1.Write()
hist_all_fCPV_1.Write()

hist_gpv.Write()
hist_gpv_1.Write()

hist_pileup.Write()
hist_pileup_1.Write()

hist_all_filt.Write()
hist_all_filt_1.Write()

hist_double_p.Write()
hist_double_1_p.Write()

outFile.Close()

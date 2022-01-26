from ROOT import TFile, TNtuple
import csv
import numpy as np

f = TFile.Open("data15.root")   #or 16,17,18

if f.IsZombie() or not f.IsOpen():
    print("Error opening file")

# Reduced variables without the jets information.
used_ntuple_variables = ["nJets_OR_TauOR", "nJets_OR_DL1r_70", "l2SS1tau", "met_met", "met_phi", \
                            "lep_Pt_0", "lep_E_0", "lep_Eta_0", "lep_Phi_0", \
                            "lep_Pt_1", "lep_E_1", "lep_Eta_1", "lep_Phi_1", \
                            "taus_pt_0", "taus_eta_0", "taus_phi_0", \
                            "HT", "taus_numTrack_0", "eventNumber","lep_ID_0","lep_isMedium_0","lep_isolationFCLoose_0","passPLIVVeryTight_0", \
                            "lep_isTightLH_0","lep_chargeIDBDTResult_recalc_rel207_tight_0","passPLIVVeryTight_0","lep_ID_1","lep_isMedium_1","lep_isolationFCLoose_1","passPLIVVeryTight_1","lep_isTightLH_1", \
                            "lep_RadiusCO_1","lep_Mtrktrk_atPV_CO_1","lep_ambiguityType_1","lep_Mtrktrk_atConvV_CO_0","lep_RadiusCO_0","lep_Mtrktrk_atPV_CO_0","lep_ambiguityType_0","lep_chargeIDBDTResult_recalc_rel207_tight_1", \
                            "dilep_type","nTaus_OR_Pt25","lep_Mtrktrk_atConvV_CO_1", \
                        ]

tree = f.nominal
tree.SetBranchStatus("*",0) 

for var_name in used_ntuple_variables:
        tree.SetBranchStatus(var_name,1)

c=0
rows = []
for event in tree:
    # The strict selection similar to the decay channel tag & num of jets requirement.
    if ((abs(event.lep_ID_0)==13 and ord(event.lep_isMedium_0) and ord(event.lep_isolationFCLoose_0) and event.passPLIVVeryTight_0) or (abs(event.lep_ID_0)==11 and ord(event.lep_isolationFCLoose_0) and ord(event.lep_isTightLH_0) and event.lep_chargeIDBDTResult_recalc_rel207_tight_0>0.7 and event.passPLIVVeryTight_0)) and ((abs(event.lep_ID_1)==13 and ord(event.lep_isMedium_1) and ord(event.lep_isolationFCLoose_1) and event.passPLIVVeryTight_1) or (abs(event.lep_ID_1)==11 and ord(event.lep_isolationFCLoose_1) and ord(event.lep_isTightLH_1) and event.lep_chargeIDBDTResult_recalc_rel207_tight_1>0.7 and event.passPLIVVeryTight_1)) and (((abs(event.lep_ID_0) == 13) or ( abs( event.lep_ID_0 ) == 11 and ord(event.lep_ambiguityType_0) == 0 and ( not ((event.lep_Mtrktrk_atPV_CO_0<0.1 and event.lep_Mtrktrk_atPV_CO_0>0) and  not (event.lep_RadiusCO_0>20 and (event.lep_Mtrktrk_atConvV_CO_0<0.1 and event.lep_Mtrktrk_atConvV_CO_0>0))) and  not (event.lep_RadiusCO_0>20 and (event.lep_Mtrktrk_atConvV_CO_0<0.1 and event.lep_Mtrktrk_atConvV_CO_0>0))))) and ((abs( event.lep_ID_1 ) == 11 and ord(event.lep_ambiguityType_1) == 0 and  not ((event.lep_Mtrktrk_atPV_CO_1<0.1 and event.lep_Mtrktrk_atPV_CO_1>0) and  not (event.lep_RadiusCO_1>20 and (event.lep_Mtrktrk_atConvV_CO_1<0.1 and event.lep_Mtrktrk_atConvV_CO_1>0))) and  not (event.lep_RadiusCO_1>20 and (event.lep_Mtrktrk_atConvV_CO_1<0.1 and event.lep_Mtrktrk_atConvV_CO_1>0))) or (abs(event.lep_ID_1) == 13))) and ord(event.nTaus_OR_Pt25)>=1 and (ord(event.nJets_OR_TauOR)>2 and ord(event.nJets_OR_DL1r_70)>0) and (event.dilep_type and event.lep_ID_0*event.lep_ID_1>0):
        c +=1
        row = []
        row.append(float(event.met_met))
        row.append(float(event.met_phi))
        row.append(float(event.lep_Pt_0))
        row.append(float(event.lep_E_0))
        row.append(float(event.lep_Eta_0))
        row.append(float(event.lep_Phi_0))
        row.append(float(event.lep_Pt_1))
        row.append(float(event.lep_E_1))
        row.append(float(event.lep_Eta_1))
        row.append(float(event.lep_Phi_1))
        row.append(float(event.taus_pt_0))
        row.append(float(event.taus_eta_0))
        row.append(float(event.taus_phi_0))
        row.append(float(event.HT))
        row.append(float(ord(event.taus_numTrack_0)))
        row.append(float(event.eventNumber))
        rows.append(row)

f = open("mc15.csv", "a")   # 16,17,18
writer = csv.writer(f)
writer.writerows(rows)
f.close()
print(c)
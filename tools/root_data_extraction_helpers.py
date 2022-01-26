from ROOT import TLorentzVector



class TruthEvent():
    def __init__(self):
        self.met_met = 0
        self.met_phi = 0
        self.usable = True
        self.main = TLorentzVector(0,0,0,0)
        self.tau_1 = TLorentzVector(0,0,0,0)
        self.tau_2 = TLorentzVector(0,0,0,0)
        self.n_1 = TLorentzVector(0,0,0,0)
        self.n_2 = TLorentzVector(0,0,0,0)
        self.hadr_tau = TLorentzVector(0,0,0,0)
        self.hadr_tau_n = TLorentzVector(0,0,0,0)
        self.hadr_tau_children = []
        self.hadr_tau_children_ids = []
        self.hadr_tau_W_q = TLorentzVector(0,0,0,0)
        self.hadr_tau_W_anti_q = TLorentzVector(0,0,0,0)
        self.lep_tau = TLorentzVector(0,0,0,0)
        self.lep_tau_n = TLorentzVector(0,0,0,0)
        self.lep_tau_children = []
        self.lep_tau_children_ids = []
        self.lep_tau_W_l = TLorentzVector(0,0,0,0)
        self.lep_tau_W_n = TLorentzVector(0,0,0,0)
        self.top = TLorentzVector(0,0,0,0)
        self.top_b = TLorentzVector(0,0,0,0)
        self.top_W = TLorentzVector(0,0,0,0)
        self.top_W_q = TLorentzVector(0,0,0,0)
        self.top_W_anti_q = TLorentzVector(0,0,0,0)
        self.top_W_l = TLorentzVector(0,0,0,0)
        self.top_W_n = TLorentzVector(0,0,0,0)
        self.anti_top = TLorentzVector(0,0,0,0)
        self.anti_top_b = TLorentzVector(0,0,0,0)
        self.anti_top_W = TLorentzVector(0,0,0,0)
        self.anti_top_W_q = TLorentzVector(0,0,0,0)
        self.anti_top_W_anti_q = TLorentzVector(0,0,0,0)
        self.anti_top_W_l = TLorentzVector(0,0,0,0)
        self.anti_top_W_n = TLorentzVector(0,0,0,0)

class RecoEvent():
    def __init__(self):
        self.hadr_tau = TLorentzVector(0,0,0,0)
        self.leps = []
        self.jets = []
        self.jets_tags = []
        self.JetsET = 0
        self.numJets25 = 0
        self.HT = 0
        self.taus_numTrack_0 = 0
        self.eventNumber = 0
        self.met_x = 0
        self.met_y = 0
        self.RunYear = 0
        self.custTrigSF_LooseID_FCLooseIso_DLT = 0 
        self.weight_pileup = 0
        self.jvtSF_customOR = 0 
        self.bTagSF_weight_DL1r_70 = 0
        self.weight_mc = 0
        self.xs = 0
        self.lep_SF_CombinedTight_0 = 0
        self.lep_SF_CombinedTight_1 = 0
        self.totalEventsWeighted= 0

    def count_numJets25_JetsET(self):
        for jet in self.jets:
            self.JetsET += jet.Et()
            if jet.Et() > 25000:
                    self.numJets25 += 1

class PermutationEvent():
    def __init__(self):
        self.top_W_q1 = TLorentzVector(0,0,0,0)
        self.top_W_q2 = TLorentzVector(0,0,0,0)
        self.top_W_q1_tag = 0
        self.top_W_q2_tag = 0
        self.top_W = TLorentzVector(0,0,0,0)
        self.top_b = TLorentzVector(0,0,0,0)
        self.anti_top_b = TLorentzVector(0,0,0,0)
        self.top_b_tag = 0
        self.anti_top_b_tag = 0
        self.main_lep = TLorentzVector(0,0,0,0)
        self.anti_top_lep = TLorentzVector(0,0,0,0)
        self.main = TLorentzVector(0,0,0,0)
        self.event_reco = RecoEvent()
        self.label = [0,0,0,0,0]
    
    def hadr_tau(self):
        return self.event_reco.hadr_tau
    
    def top(self):
        return self.top_W + self.top_b
    
    def anti_top(self):
        return self.anti_top_lep + self.anti_top_b
    
    def visible_main(self):
        return self.main_lep + self.event_reco.hadr_tau

    def get_reco_jets(self):
        return [self.top_b, self.anti_top_b, self.top_W_q1, self.top_W_q2]
        #return [self.top_b, self.anti_top_b, self.top_W]
    
    def get_reco_leps(self):
        return [self.main_lep, self.anti_top_lep]

    def get_tags(self):
        return [self.top_b_tag, self.anti_top_b_tag, self.top_W_q1_tag, self.top_W_q2_tag]

    def get_delta_rs(self):
        return[dist(self.top_W_q1, self.top_W_q2), \
               dist(self.hadr_tau(), self.main_lep), \
               dist(self.hadr_tau(), self.anti_top_lep), \
               dist(self.anti_top_b, self.main_lep), \
               dist(self.anti_top_b, self.anti_top_lep), \
               dist(self.top_b, self.top_W), \
               dist(self.top(), self.anti_top())]

class OutputEvent(PermutationEvent):
    def __init__(self):
        super().__init__()
        
def dist(p1,p2):
    return p1.DeltaR(p2)

def generate_jets_mask(event):
    jet_pt0 = event.jet_tauOR_pt0
    jet_pt1 = event.jet_tauOR_pt1
    jet_pt2 = event.jet_tauOR_pt2
    jet_pt3 = event.jet_tauOR_pt3
    jet_pt4 = event.jet_tauOR_pt4
    jet_pt5 = event.jet_tauOR_pt5
    jet_pt6 = event.jet_tauOR_pt6
    jet_pt7 = event.jet_tauOR_pt7
    jets_pts = [jet_pt0 != -99.0,jet_pt1 != -99.0,jet_pt2 != -99.0,jet_pt3 != -99.0,jet_pt4 != -99.0,jet_pt5 != -99.0,jet_pt6 != -99.0,jet_pt7 != -99.0]
    return jets_pts

def get_jets(event,jets_mask):
    q1 = TLorentzVector(0,0,0,0)    # q1 and q2 chosen in order
    q2 = TLorentzVector(0,0,0,0)
    q3 = TLorentzVector(0,0,0,0)
    q4 = TLorentzVector(0,0,0,0)
    q5 = TLorentzVector(0,0,0,0)
    q6 = TLorentzVector(0,0,0,0)
    q7 = TLorentzVector(0,0,0,0)
    q8 = TLorentzVector(0,0,0,0)

    indices = [i for i, x in enumerate(jets_mask) if x]

    if len(indices) > 0:
        jet_pt0 = getattr(event, "jet_tauOR_pt"+str(indices[0]))
        jet_eta0 = getattr(event, "jet_tauOR_eta"+str(indices[0]))
        jet_phi0 = getattr(event, "jet_tauOR_phi"+str(indices[0]))
        jet_E0 = getattr(event, "jet_tauOR_E"+str(indices[0]))
        q1.SetPtEtaPhiE(jet_pt0,jet_eta0,jet_phi0,jet_E0)
    if len(indices) > 1:
        jet_pt1 = getattr(event, "jet_tauOR_pt"+str(indices[1]))
        jet_eta1 = getattr(event, "jet_tauOR_eta"+str(indices[1]))
        jet_phi1 = getattr(event, "jet_tauOR_phi"+str(indices[1]))
        jet_E1 = getattr(event, "jet_tauOR_E"+str(indices[1]))
        q2.SetPtEtaPhiE(jet_pt1,jet_eta1,jet_phi1,jet_E1)
    if len(indices) > 2:
        jet_pt2 = getattr(event, "jet_tauOR_pt"+str(indices[2]))
        jet_eta2 = getattr(event, "jet_tauOR_eta"+str(indices[2]))
        jet_phi2 = getattr(event, "jet_tauOR_phi"+str(indices[2]))
        jet_E2 = getattr(event, "jet_tauOR_E"+str(indices[2]))
        q3.SetPtEtaPhiE(jet_pt2,jet_eta2,jet_phi2,jet_E2)
    if len(indices) > 3:
        jet_pt3 = getattr(event, "jet_tauOR_pt"+str(indices[3]))
        jet_eta3 = getattr(event, "jet_tauOR_eta"+str(indices[3]))
        jet_phi3 = getattr(event, "jet_tauOR_phi"+str(indices[3]))
        jet_E3 = getattr(event, "jet_tauOR_E"+str(indices[3]))
        q4.SetPtEtaPhiE(jet_pt3,jet_eta3,jet_phi3,jet_E3)
    if len(indices) > 4:
        jet_pt4 = getattr(event, "jet_tauOR_pt"+str(indices[4]))
        jet_eta4 = getattr(event, "jet_tauOR_eta"+str(indices[4]))
        jet_phi4 = getattr(event, "jet_tauOR_phi"+str(indices[4]))
        jet_E4 = getattr(event, "jet_tauOR_E"+str(indices[4]))
        q5.SetPtEtaPhiE(jet_pt4,jet_eta4,jet_phi4,jet_E4)
    if len(indices) > 5:
        jet_pt5 = getattr(event, "jet_tauOR_pt"+str(indices[5]))
        jet_eta5 = getattr(event, "jet_tauOR_eta"+str(indices[5]))
        jet_phi5 = getattr(event, "jet_tauOR_phi"+str(indices[5]))
        jet_E5 = getattr(event, "jet_tauOR_E"+str(indices[5]))
        q6.SetPtEtaPhiE(jet_pt5,jet_eta5,jet_phi5,jet_E5)
    if len(indices) > 6:
        jet_pt6 = getattr(event, "jet_tauOR_pt"+str(indices[6]))
        jet_eta6 = getattr(event, "jet_tauOR_eta"+str(indices[6]))
        jet_phi6 = getattr(event, "jet_tauOR_phi"+str(indices[6]))
        jet_E6 = getattr(event, "jet_tauOR_E"+str(indices[6]))
        q7.SetPtEtaPhiE(jet_pt6,jet_eta6,jet_phi6,jet_E6)
    if len(indices) > 7:
        jet_pt7 = getattr(event, "jet_tauOR_pt"+str(indices[7]))
        jet_eta7 = getattr(event, "jet_tauOR_eta"+str(indices[7])) 
        jet_phi7 = getattr(event, "jet_tauOR_phi"+str(indices[7]))   
        jet_E7 = getattr(event, "jet_tauOR_E"+str(indices[7]))     
        q8.SetPtEtaPhiE(jet_pt7,jet_eta7,jet_phi7,jet_E7)
    
    fake_jet = TLorentzVector()
    fake_jet.SetPtEtaPhiE(0,0,0,0)

    if len(indices) == 0:
        return([fake_jet,fake_jet,fake_jet,fake_jet])
    if len(indices) == 1:
        return([q1,fake_jet,fake_jet,fake_jet])
    if len(indices) == 2:
        return([q1,q2,fake_jet,fake_jet])
    if len(indices) == 3:
        return([q1,q2,q3,fake_jet])
    if len(indices) == 4:
        return([q1,q2,q3,q4])
    if len(indices) == 5:
        return([q1,q2,q3,q4,q5])
    if len(indices) == 6:
        return([q1,q2,q3,q4,q5,q6])
    if len(indices) == 7:
        return([q1,q2,q3,q4,q5,q6,q7])
    if len(indices) == 8:
        return([q1,q2,q3,q4,q5,q6,q7,q8])

def get_btags(event,jets):
    btags_ret = []
    btags = [event.jet_pseudoscore_DL1r0, event.jet_pseudoscore_DL1r1, event.jet_pseudoscore_DL1r2, event.jet_pseudoscore_DL1r3, event.jet_pseudoscore_DL1r4]
    btag_index_delta = 0 

    for i in range(len(jets)):
        jet = jets[i]

        btag = 1
        if i-btag_index_delta <= 4:
            btag = btags[i-btag_index_delta]

            if btag < 0:
                btag = 1

        if jet.Pt() != 0.0:
            btags_ret.append(btag)
        else:
            btag_index_delta += 1
            btags_ret.append(0)

    return btags_ret
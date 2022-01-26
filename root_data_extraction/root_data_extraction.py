from tools.root_data_extraction_helpers import *
from itertools import permutations, combinations
from ROOT import TLorentzVector, TFile
import numpy as np
import csv
from tools.constants import *
from os import listdir
from os.path import isfile,isdir,join

def extract_main_particle_from_truth_event(event, main_particle_id):

    truth_information = TruthEvent()

    if main_particle_id == -1:    # Means there is no main particle.
        truth_information.main = TLorentzVector(0,0,0,0)
        return truth_information

    elif main_particle_id == 24:    # W boson.

        """ The truth information is contained in vectors, where each particle is represented
        by an index with which we access each array to obtain the particle's information. """
        """ Transform vector of vectors to a single dimension array. """
        parent_arrays = event.m_truth_parents
        parents = []
        for parent_array in parent_arrays:
            try:
                p = parent_array[0]
            except:
                p = None
            parents.append(p)

        """ BARCODE = unique identificator
            pdgId   = particle type identificator (e.g. pdgId = 24 means particle is W boson) """
        barcodes = event.m_truth_barcode
        ids = event.m_truth_pdgId

        pts = event.m_truth_pt
        etas = event.m_truth_eta
        phis = event.m_truth_phi
        es = event.m_truth_e

        for index,id in enumerate(ids):
            if abs(id) == 24 and parents[index] not in barcodes:      # W boson with no parent is the one we are looking for.
                main_W_boson = TLorentzVector()
                main_W_boson.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                truth_information.main = main_W_boson
                return truth_information

    else:
        """ Higgs or Z boson. Similar to W boson, but without the no-parent clause. """
        parent_arrays = event.m_truth_parents
        parents = []
        for parent_array in parent_arrays:
            try:
                p = parent_array[0]
            except:
                p = None
            parents.append(p)

        barcodes = event.m_truth_barcode
        ids = event.m_truth_pdgId

        pts = event.m_truth_pt
        etas = event.m_truth_eta
        phis = event.m_truth_phi
        es = event.m_truth_e

        for index,id in enumerate(ids):
            if abs(id) == main_particle_id:
                main_H_or_Z_boson = TLorentzVector()
                main_H_or_Z_boson.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                truth_information.main = main_H_or_Z_boson
                return truth_information

def extract_from_truth_event(event, main_particle_id):

    """ Helper variables for particles with children. """
    main_barcode = 0
    tau_1_barcode = 0
    tau_2_barcode = 0
    hadr_tau_n_barcode = 0
    lep_tau_n_barcode = 0
    top_barcode = 0
    top_W_barcode = 0
    top_W_tau_barcode = 0
    anti_top_barcode = 0
    anti_top_W_barcode = 0
    anti_top_W_tau_barcode = 0

    """ Container object. """
    truth_information = TruthEvent()
    
    """ The truth information is contained in vectors, where each particle is represented
        by an index with which we access each array to obtain the particle's information. """
    """ Transform vector of vectors to a single dimension array. """
    parent_arrays = event.m_truth_parents
    parents = []
    for parent_array in parent_arrays:
        try:
            p = parent_array[0]
        except:
            p = None
        parents.append(p)

    """ BARCODE = unique identificator
        pdgId   = particle type identificator (e.g. pdgId = 25 means particle is Higgs boson) """

    barcodes = event.m_truth_barcode
    ids = event.m_truth_pdgId

    pts = event.m_truth_pt
    etas = event.m_truth_eta
    phis = event.m_truth_phi
    es = event.m_truth_e

    solved_barcodes = set()

    """ Iterate through the particles in direction from origin particles to final state particles. """

    """ First iteration finds main particle (e.g. Higgs) and two top quarks. """
    for index,id in enumerate(ids):
        if abs(id) == main_particle_id:
            main_barcode = barcodes[index]
            truth_information.main.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
        elif id == 6:
            top_barcode = barcodes[index]
            truth_information.top.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
        elif id == -6:
            anti_top_barcode = barcodes[index]
            truth_information.anti_top.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])

    """ Second iteration finds top Ws, top b quarks and Higgs taus. """
    for index,parent_and_id in enumerate(list(zip(parents,ids))):
        parent = parent_and_id[0]
        id = parent_and_id[1]

        if parent == main_barcode: # Taus can't actually be distinguished at this point (hadr or lep decay). We use the generic tau1 and tau2 members.
            if id == 15:
                tau_1_barcode = barcodes[index]
                truth_information.tau_1.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
            elif id == -15:
                tau_2_barcode = barcodes[index]
                truth_information.tau_2.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])

        elif parent == top_barcode:
            if id == 5:
                top_b_barcode = barcodes[index]
                truth_information.top_b.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
            elif id == 24:
                top_W_barcode = barcodes[index]
                truth_information.top_W.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])

        elif parent == anti_top_barcode:
            if id == -5:
                anti_top_b_barcode = barcodes[index]
                truth_information.anti_top_b.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
            elif id == -24:
                anti_top_W_barcode = barcodes[index]
                truth_information.anti_top_W.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
    
    """ The third and last iteration is more complicated due to the fact, that the Ws and taus in the data can decay
        into another tau or W (even multiple times) before decaying into their children particles. We can even have
        We have to iterate through these multiple-layered Ws and taus until we get to the children. """
    phase_solved = False

    while not phase_solved:
        phase_solved = True

        for index,parent_and_id in enumerate(list(zip(parents,ids))):
            parent = parent_and_id[0]
            id = parent_and_id[1]

            if parent == top_W_barcode:
                if id == 1 or id == 2 or id == 3 or id == 4:
                    truth_information.top_W_q.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif id == -1 or id == -2 or id == -3 or id == -4:
                    truth_information.top_W_anti_q.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif abs(id) == 11 or abs(id) == 13:
                    truth_information.top_W_l.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif abs(id) == 12 or abs(id) == 14:
                    truth_information.top_W_n.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif abs(id) == 15:
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        phase_solved = False
                        top_W_tau_barcode = barcodes[index]
                elif abs(id) == 16:
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        neutrino = TLorentzVector()
                        neutrino.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                        truth_information.top_W_n += neutrino 
                elif id == 24:
                    phase_solved = False
                    top_W_barcode = barcodes[index]
                    truth_information.top_W.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                else:
                    truth_information.usable = False

            elif parent == tau_1_barcode:
                if abs(id) == 111 or abs(id) == 211 or abs(id) == 321 or abs(id) == 311 \
                or abs(id) == 221 or abs(id) == 310 or abs(id) == 223 or abs(id) == 323 \
                or abs(id) == 130 or abs(id) == 1 or abs(id) == 2 or abs(id) == 3 or abs(id) == 4 :
                    truth_information.hadr_tau = truth_information.tau_1
                    truth_information.hadr_tau_n = truth_information.n_1
                    child = TLorentzVector(0,0,0,0)
                    child.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        truth_information.hadr_tau_children.append(child)
                        truth_information.hadr_tau_children_ids.append(id)
                elif abs(id) == 11 or abs(id) == 12 or abs(id) == 13 or abs(id) == 14:
                    truth_information.lep_tau = truth_information.tau_1
                    truth_information.lep_tau_n = truth_information.n_1
                    child = TLorentzVector(0,0,0,0)
                    child.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        if abs(id) == 11 or abs(id) == 13:
                            truth_information.lep_tau_W_l = child
                        elif abs(id) == 12 or abs(id) == 14:
                            truth_information.lep_tau_W_n = child
                elif abs(id) == 16:
                    truth_information.n_1.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif abs(id) == 15:
                    phase_solved = False
                    tau_1_barcode = barcodes[index]
                    truth_information.tau_1.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])

            elif parent == tau_2_barcode:
                if abs(id) == 111 or abs(id) == 211 or abs(id) == 321 or abs(id) == 311 \
                or abs(id) == 221 or abs(id) == 310 or abs(id) == 223 or abs(id) == 323 \
                or abs(id) == 130 or abs(id) == 1 or abs(id) == 2 or abs(id) == 3 or abs(id) == 4 :
                    truth_information.hadr_tau = truth_information.tau_2
                    truth_information.hadr_tau_n = truth_information.n_2
                    child = TLorentzVector(0,0,0,0)
                    child.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        truth_information.hadr_tau_children.append(child)
                        truth_information.hadr_tau_children_ids.append(id)
                elif abs(id) == 11 or abs(id) == 12 or abs(id) == 13 or abs(id) == 14:
                    truth_information.lep_tau = truth_information.tau_2
                    truth_information.lep_tau_n = truth_information.n_2
                    child = TLorentzVector(0,0,0,0)
                    child.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        if abs(id) == 11 or abs(id) == 13:
                            truth_information.lep_tau_W_l = child
                        elif abs(id) == 12 or abs(id) == 14:
                            truth_information.lep_tau_W_n = child
                elif abs(id) == 16:
                    truth_information.n_2.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif abs(id) == 15:
                    phase_solved = False
                    tau_2_barcode = barcodes[index]
                    truth_information.tau_2.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])

            elif parent == anti_top_W_barcode:
                if id == 1 or id == 2 or id == 3 or id == 4:
                    truth_information.anti_top_W_q.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif id == -1 or id == -2 or id == -3 or id == -4:
                    truth_information.anti_top_W_anti_q.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif abs(id) == 11 or abs(id) == 13:
                    truth_information.anti_top_W_l.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif abs(id) == 12 or abs(id) == 14:
                    truth_information.anti_top_W_n.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif abs(id) == 15:
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        phase_solved = False
                        anti_top_W_tau_barcode = barcodes[index]
                elif abs(id) == 16:
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        neutrino = TLorentzVector()
                        neutrino.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                        truth_information.anti_top_W_n += neutrino
                elif id == -24:
                    phase_solved = False
                    anti_top_W_barcode = barcodes[index]
                    truth_information.anti_top_W.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                else:
                    truth_information.usable = False
            
            elif parent == top_W_tau_barcode:  # If this tau decays hadronically we cannot use the event for particle assignment.
                if abs(id) == 11 or abs(id) == 13:
                    truth_information.top_W_l.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif abs(id) == 12 or abs(id) == 14:
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        neutrino = TLorentzVector()
                        neutrino.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                        truth_information.top_W_n += neutrino
                elif abs(id) == 16:
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        neutrino = TLorentzVector()
                        neutrino.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                        truth_information.top_W_n += neutrino
                elif abs(id) == 15:
                    phase_solved = False
                    top_W_tau_barcode = barcodes[index]

            elif parent == anti_top_W_tau_barcode:  # If this tau decays hadronically we cannot use the event for particle assignment.
                if abs(id) == 11 or abs(id) == 13:
                    truth_information.anti_top_W_l.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                elif abs(id) == 12 or abs(id) == 14:
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        neutrino = TLorentzVector()
                        neutrino.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                        truth_information.anti_top_W_n += neutrino
                elif abs(id) == 16:
                    if not barcodes[index] in solved_barcodes:
                        solved_barcodes.add(barcodes[index])
                        neutrino = TLorentzVector()
                        neutrino.SetPtEtaPhiE(pts[index],etas[index],phis[index],es[index])
                        truth_information.anti_top_W_n += neutrino
                elif abs(id) == 15:
                    phase_solved = False
                    anti_top_W_tau_barcode = barcodes[index]

    return truth_information

def extract_from_reco_event(event):
    reco_information = RecoEvent()

    """ Jets have to be extracted with the use of a helper mask because of the nature of the data.
        In lated productions of the root ntuples the jets are stored differently
        and their extraction would have to be changed. """
    jets_mask = generate_jets_mask(event)
    reco_information.jets = get_jets(event,jets_mask)
    reco_information.jets_tags = get_btags(event,reco_information.jets)
    reco_information.count_numJets25_JetsET()

    reco_information.taus_numTrack_0 = ord(event.taus_numTrack_0)*10
    reco_information.HT = event.HT
    reco_information.eventNumber = event.eventNumber

    """ Almost always there are two detected leptons but in very rare cases there is also a third jet. """
    l1 = TLorentzVector()
    l2 = TLorentzVector()
    l1.SetPtEtaPhiE(event.lep_Pt_0, event.lep_Eta_0, event.lep_Phi_0, event.lep_E_0)
    l2.SetPtEtaPhiE(event.lep_Pt_1, event.lep_Eta_1, event.lep_Phi_1, event.lep_E_1)
    
    if not np.isclose(event.lep_Pt_2, 0.0):
        l3 = TLorentzVector()
        l3.SetPtEtaPhiE(event.lep_Pt_2, event.lep_Eta_2, event.lep_Phi_2, event.lep_E_2)
        reco_information.leps = [l1,l2,l3]
    else:
        reco_information.leps = [l1,l2]

    hadr_tau = TLorentzVector()
    hadr_tau.SetPtEtaPhiE(event.taus_pt_0, event.taus_eta_0, event.taus_phi_0, event.taus_E_0)
    reco_information.hadr_tau = hadr_tau

    """ Transform missing transverse energy from transverse energy and phi representation to transverse momentum representation. """
    met_x = event.met_met * np.cos(event.met_phi)
    met_y = event.met_met * np.sin(event.met_phi)
    reco_information.met_x = met_x
    reco_information.met_y = met_y

    """ Add variables for sample weights. """
    """ reco_information.RunYear = event.RunYear
    reco_information.custTrigSF_LooseID_FCLooseIso_DLT = event.custTrigSF_LooseID_FCLooseIso_DLT
    reco_information.weight_pileup = event.weight_pileup
    reco_information.jvtSF_customOR = event.jvtSF_customOR
    reco_information.bTagSF_weight_DL1r_70 = event.bTagSF_weight_DL1r_70
    reco_information.weight_mc = event.weight_mc
    reco_information.xs = event.xs
    reco_information.lep_SF_CombinedTight_0 = event.lep_SF_CombinedTight_0
    reco_information.lep_SF_CombinedTight_1 = event.lep_SF_CombinedTight_1
    reco_information.totalEventsWeighted = event.totalEventsWeighted """

    return reco_information
    #reco_data.append([jets, leps, [hadr_tau], [met_x, met_y, JestET, numJets25, taus_numTrack_0, HT],[event.eventNumber]]) # check indices if changed !!!!!!!

def get_label(reco_jets,reco_leps,true_jets,true_leps):
        label = []

        for i in range(2):
            if not dist(reco_jets[i],true_jets[i]) < JETS_THRESHOLD or reco_jets[i].Px() == 0:
                label.append(0)
            else:
                label.append(1)

        if reco_jets[2].Px() == 0 or reco_jets[3].Px() == 0:
            label.append(0)
        else:
            if dist(reco_jets[2],true_jets[2]) < JETS_THRESHOLD and dist(reco_jets[3],true_jets[3]) < JETS_THRESHOLD:
                label.append(1)
            elif dist(reco_jets[2],true_jets[3]) < JETS_THRESHOLD and dist(reco_jets[3],true_jets[2]) < JETS_THRESHOLD:
                label.append(1)
            else:
                label.append(0)

        for j in range(len(reco_leps)):
            if not dist(reco_leps[j],true_leps[j]) < LEPTONS_THRESHOLD:
                label.append(0)
            else:
                label.append(1)

        return label

def generate_permutations(event_reco, event_truth):

    """ For an event create all possible permutations of jets and leptons at positions
        top b jet, anti-top b-jet, top W for jets and main (Higgs) lepton, anti-top lepton
        for leptons. """
    event_permutations = []

      # True jets and leptons are needed for creating labels for particle assignment later in this func.
    if event_truth.anti_top_W_l.Mag() != 0:                                     # Whichever top decays leptonically.
        true_leps = [event_truth.lep_tau_W_l, event_truth.anti_top_W_l]    
        true_jets = [event_truth.top_b, event_truth.anti_top_b, event_truth.top_W_q, event_truth.top_W_anti_q]    
    else:
        true_leps = [event_truth.lep_tau_W_l, event_truth.top_W_l]         
        true_jets = [event_truth.anti_top_b, event_truth.top_b, event_truth.anti_top_W_q, event_truth.anti_top_W_anti_q]

    """ First all possible pair combinations of jets. These pairs will make up the top W boson. """
    W_pairs = list(combinations(zip(event_reco.jets,event_reco.jets_tags), r=2))                  # All possible Ws.
    W_pairs_pxs = [[pair[0][0].Px(),pair[1][0].Px()] for pair in W_pairs]     # Distinguish between jets based on their Px.

    rest_of_jets = []
    for pxs in W_pairs_pxs:
        rest_of_jets.append([jet for jet in zip(event_reco.jets,event_reco.jets_tags) if jet[0].Px() not in pxs]) # All jets not already in the W pair
        while len(rest_of_jets[-1]) < 2:
            rest_of_jets[-1] = rest_of_jets[-1] + [(TLorentzVector(0,0,0,0),0)]             # If less than two jets, supplement with "zero jets".

    """ Second all permutations of the rest of the jets not in the W."""
    for j in range(len(W_pairs)):
        W_pair = W_pairs[j]
        rest = rest_of_jets[j]
        W_summed = (W_pair[0][0]+W_pair[1][0])

        jets_permutations = list(permutations(rest, r=2))
        
        """ Lastly all lepton permutations."""
        for jet_perm in jets_permutations:
            jet_perm = list(jet_perm)

            leps_permutations = list(permutations(event_reco.leps, r=2))

            for lep_perm in leps_permutations:
                lep_perm = list(lep_perm)

                output_event = PermutationEvent()    # Container object.

                output_event.top_W_q1 = W_pair[0][0] if W_pair[0][0].E() > W_pair[1][0].E() else W_pair[1][0]
                output_event.top_W_q2 = W_pair[0][0] if W_pair[0][0].E() <= W_pair[1][0].E() else W_pair[1][0]
                output_event.top_W_q1_tag = W_pair[0][1] if W_pair[0][0].E() > W_pair[1][0].E() else W_pair[1][1]
                output_event.top_W_q2_tag = W_pair[0][1] if W_pair[0][0].E() <= W_pair[1][0].E() else W_pair[1][1]

                output_event.top_W = W_summed

                output_event.top_b = jet_perm[0][0]
                output_event.anti_top_b = jet_perm[1][0]
                output_event.top_b_tag = jet_perm[0][1]
                output_event.anti_top_b_tag = jet_perm[1][1]

                output_event.main_lep = lep_perm[0]
                output_event.anti_top_lep = lep_perm[1]

                output_event.main = event_truth.main
                output_event.nu1 = event_truth.main
                output_event.nu2 = event_truth.main
                output_event.nu3 = event_truth.main
                output_event.nu4 = event_truth.main

                output_event.event_reco = event_reco

                
                label = get_label(output_event.get_reco_jets(), output_event.get_reco_leps(), true_jets, true_leps)
                output_event.label = label

                event_permutations.append(output_event)

    return event_permutations

def generate_variables_to_write(event_permutations):
    events_to_write = []
    labels_to_write = []

    """ First fetch all the variables that will be constant for each permutation. """
    hadr_tau = event_permutations[0].hadr_tau()
    main = event_permutations[0].main
    eventNumber = event_permutations[0].event_reco.eventNumber
    JetsET = event_permutations[0].event_reco.JetsET
    numJets25 = event_permutations[0].event_reco.numJets25
    HT = event_permutations[0].event_reco.HT
    taus_numTrack_0 = event_permutations[0].event_reco.taus_numTrack_0
    met_x = event_permutations[0].event_reco.met_x
    met_y = event_permutations[0].event_reco.met_y

    num_permutations = len(event_permutations)

    """ For each permutation get all the variables we want.
        Unfortunately cannot be done with vectors, because we are using the ROOT TLorentzVector objects
        and other ROOT methods.
        But the runtime is fine for a one-time operation... """
    for i in range(num_permutations):
        event_permutation = event_permutations[i]
        event_to_write = [eventNumber]                   # event number serves as event ID
        event_to_write += [event_permutation.top_b.Px(), event_permutation.top_b.Py(), event_permutation.top_b.Pz(), event_permutation.top_b.E()] # top b
        event_to_write += [event_permutation.anti_top_b.Px(), event_permutation.anti_top_b.Py(), event_permutation.anti_top_b.Pz(), event_permutation.anti_top_b.E()] # antitop b
        event_to_write += [event_permutation.top_W.Px(), event_permutation.top_W.Py(), event_permutation.top_W.Pz(), event_permutation.top_W.E()] # W
        event_to_write += [event_permutation.main_lep.Px(), event_permutation.main_lep.Py(), event_permutation.main_lep.Pz(), event_permutation.main_lep.E()] # main lepton
        event_to_write += [event_permutation.anti_top_lep.Px(), event_permutation.anti_top_lep.Py(), event_permutation.anti_top_lep.Pz(), event_permutation.anti_top_lep.E()] # antitop lepton
        event_to_write += [event_permutation.top_W_q1.Px(), event_permutation.top_W_q1.Py(), event_permutation.top_W_q1.Pz(), event_permutation.top_W_q1.E()] # q1
        event_to_write += [event_permutation.top_W_q2.Px(), event_permutation.top_W_q2.Py(), event_permutation.top_W_q2.Pz(), event_permutation.top_W_q2.E()] # q2
        event_to_write += [hadr_tau.Px(), hadr_tau.Py(), hadr_tau.Pz(), hadr_tau.E()] # tau
        event_to_write += [event_permutation.top().Px(), event_permutation.top().Py(), event_permutation.top().Pz(), event_permutation.top().E()]
        event_to_write += [event_permutation.anti_top().Px(), event_permutation.anti_top().Py(), event_permutation.anti_top().Pz(), event_permutation.anti_top().E()]
        event_to_write += [event_permutation.visible_main().Px(), event_permutation.visible_main().Py(), event_permutation.visible_main().Pz(), event_permutation.visible_main().E()]
        event_to_write += event_permutation.get_tags()
        event_to_write += [event_permutation.top_b.Mag(), event_permutation.anti_top_b.Mag(), event_permutation.top_W.Mag(), \
                           event_permutation.main_lep.Mag(), event_permutation.anti_top_lep.Mag(), event_permutation.top_W_q1.Mag(), \
                           event_permutation.top_W_q2.Mag(), event_permutation.top().Mag(), \
                           event_permutation.anti_top().Mag(),event_permutation.visible_main().Mag()]
        event_to_write += event_permutation.get_delta_rs()
        event_to_write += [met_x, met_y, JetsET, numJets25, HT, taus_numTrack_0, num_permutations]
        event_to_write += [main.Px(), main.Py(), main.Pz(), main.E()]   # higgs
        events_to_write.append(event_to_write)
        labels_to_write.append(event_permutation.label)
    return events_to_write, labels_to_write

def write_events_to_file(id,ids,labels_dict,labels_to_write,data_to_write,output_folder):
    id = id+1
    ids.append(id)
    labels_dict[str(id)] = (np.array(labels_to_write).flatten()).tolist()
    np.save(output_folder+"/data/"+str(id), data_to_write)
    return id,ids,labels_dict

def extract_data(file_path, file_id, output_folder, distribution, particle_assignment):

    if distribution == "ttH":
        main_particle_id = 25
    elif distribution == "ttZ":
        main_particle_id = 23
    elif distribution == "ttW":
        main_particle_id = 24
    elif distribution == "tt":
        main_particle_id = -1

    f = TFile.Open(file_path)

    if f.IsZombie() or not f.IsOpen():
        print("Error opening file")

    tree = f.nominal

    used_ntuple_variables = ["jet_pseudoscore_DL1r0","jet_pseudoscore_DL1r1","jet_pseudoscore_DL1r2","jet_pseudoscore_DL1r3","jet_pseudoscore_DL1r4", \
                            "nJets_OR_TauOR", "nJets_OR_DL1r_70", "l2SS1tau", "met_met", "met_phi", \
                            "lep_Pt_0", "lep_E_0", "lep_Eta_0", "lep_Phi_0", \
                            "lep_Pt_1", "lep_E_1", "lep_Eta_1", "lep_Phi_1", \
                            "lep_Pt_2", "lep_E_2", "lep_Eta_2", "lep_Phi_2", \
                            "taus_pt_0", "taus_E_0", "taus_eta_0", "taus_phi_0", \
                            "jet_tauOR_pt0", "jet_tauOR_E0", "jet_tauOR_eta0", "jet_tauOR_phi0", \
                            "jet_tauOR_pt1", "jet_tauOR_E1", "jet_tauOR_eta1", "jet_tauOR_phi1", \
                            "jet_tauOR_pt2", "jet_tauOR_E2", "jet_tauOR_eta2", "jet_tauOR_phi2", \
                            "jet_tauOR_pt3", "jet_tauOR_E3", "jet_tauOR_eta3", "jet_tauOR_phi3", \
                            "jet_tauOR_pt4", "jet_tauOR_E4", "jet_tauOR_eta4", "jet_tauOR_phi4", \
                            "jet_tauOR_pt5", "jet_tauOR_E5", "jet_tauOR_eta5", "jet_tauOR_phi5", \
                            "jet_tauOR_pt6", "jet_tauOR_E6", "jet_tauOR_eta6", "jet_tauOR_phi6", \
                            "jet_tauOR_pt7", "jet_tauOR_E7", "jet_tauOR_eta7", "jet_tauOR_phi7", \
                            "m_truth_m", "m_truth_pt", "m_truth_eta", "m_truth_phi", "m_truth_e", \
                            "m_truth_pdgId", "m_truth_status", "m_truth_barcode", "m_truth_children", "m_truth_parents", \
                            "HT", "taus_numTrack_0", "eventNumber", \
                            "RunYear", "custTrigSF_LooseID_FCLooseIso_DLT", "weight_pileup", "jvtSF_customOR", \
                            "bTagSF_weight_DL1r_70", "weight_mc", "xs", "lep_SF_CombinedTight_0", "lep_SF_CombinedTight_1"]

    """ Leave all unused variables with status 0 to speed up iterating through tree. """
    tree.SetBranchStatus("*",0) 
    for var_name in used_ntuple_variables:
        tree.SetBranchStatus(var_name,1)

    """ Iterate through all l2SS1tauhad events in the tree and return truth and detector simulation (reco)
        information in the form of an array of container objects. """
    
    events_truth = []
    events_reco = []
    counter = 0
    for event in tree:
        counter += 1
        if counter%1000 == 0:
            print(counter)  # Just for orientation.
        if event.l2SS1tau and ord(event.nJets_OR_TauOR) > 2 and ord(event.nJets_OR_DL1r_70) > 0:  # Only use events with positive decay channel tag & number of jets selection.
            if particle_assignment:
                events_truth.append(extract_from_truth_event(event, main_particle_id))  # For particle assignment training we need complete true particle tree.
            else:
                events_truth.append(extract_main_particle_from_truth_event(event, main_particle_id))    # For mass reco we need only true main particle.
            events_reco.append(extract_from_reco_event(event))

    print(len(events_reco))

    """ Check if all the necessary truth variables are present in extracted events
        (only applies to particle assignment). If an event is misidentified with the l2SS1tauhad tag,
        it can be missing some truth particles we need for the particle assignment
        and thus it cannot be used in it (for training). We delete all such events. """

    if particle_assignment:
        deleted = 0
        events_truth_filtered = []

        for i in range(len(events_truth)):
            event_truth = events_truth[i]
            unusable = False

            if  event_truth.top_b.Mag() == 0 or \
                event_truth.anti_top_b.Mag() == 0 or \
                event_truth.lep_tau_W_l.Mag() == 0:
                unusable = True

            elif event_truth.anti_top_W_l.Mag() == 0 and \
                event_truth.top_W_l.Mag() == 0:
                unusable = True

            elif event_truth.top_W.Mag() == 0 and \
                event_truth.anti_top_W.Mag() == 0:
                unusable = True
            
            elif (event_truth.top_W_anti_q.Mag() == 0 or event_truth.top_W_q.Mag() == 0) and \
                (event_truth.anti_top_W_anti_q.Mag() == 0 or event_truth.anti_top_W_q.Mag() == 0):
                unusable = True

            if unusable:
                index = i-deleted
                events_reco.pop(index)
                deleted += 1

            else:
                events_truth_filtered.append(event_truth)

        events_truth = events_truth_filtered

    else:
        deleted = 0
        events_truth_filtered = []

        for i in range(len(events_truth)):
            event_truth = events_truth[i]
            unusable = False

            if event_truth is None:
                unusable = True

            if unusable:
                index = i-deleted
                events_reco.pop(index)
                deleted += 1

            else:
                events_truth_filtered.append(event_truth)

        events_truth = events_truth_filtered

    print(len(events_reco))

    """ For each event generate all permutations of jets and leptons, which is necessary for particle assignment.
        Event if the production has particle_assignment False, we will want to have its data in the same format,
        so we do not worry about that here and we will ignore its labels later on. 
        
        Then generate all the other variables that we want to use as input for the neural network.
        This includes variables such as mass of a particle or delta R between two selected particles etc.
        
        Finally write to output files."""
    
    labels_dict = {}
    train_ids = []
    val_ids = []
    test_ids = []
    train_data_to_write = []
    val_data_to_write = []
    test_data_to_write = []
    train_labels_to_write = []
    val_labels_to_write = []
    test_labels_to_write = []

    id = int(1e06 * file_id + 1e04 * (main_particle_id if main_particle_id>0 else 42))
    
    for i in range(len(events_reco)):
        event_reco = events_reco[i]
        event_truth = events_truth[i]
        event_permutations = generate_permutations(event_reco, event_truth)
        events_to_write, labels_to_write = generate_variables_to_write(event_permutations)

        if i%10 == 9:       # Ten percent of events will be used for validation.
            val_data_to_write += events_to_write
            val_labels_to_write += labels_to_write

        elif i%10 == 4:     # Ten percent of events will be used for testing.
            test_data_to_write += events_to_write
            test_labels_to_write += labels_to_write

        else:               # Eighty percent of events will be used for training.
            train_data_to_write += events_to_write
            train_labels_to_write += labels_to_write
        
        """ We separate the data into files of approx. 1680 permutations. """
        if len(train_data_to_write) >= 1680:    
            id, train_ids, labels_dict = write_events_to_file(id, train_ids, labels_dict, train_labels_to_write, train_data_to_write, output_folder)

            train_labels_to_write = []
            train_data_to_write = []
        
        if len(val_data_to_write) >= 1680:
            id, val_ids, labels_dict = write_events_to_file(id, val_ids, labels_dict, val_labels_to_write, val_data_to_write, output_folder)

            val_labels_to_write = []
            val_data_to_write = []
        
        if len(test_data_to_write) >= 1680:
            id, test_ids, labels_dict = write_events_to_file(id, test_ids, labels_dict, test_labels_to_write, test_data_to_write, output_folder)
            
            test_labels_to_write = []
            test_data_to_write = []

    """ After all events have been processed, write any remains into files. """
    if len(train_data_to_write) != 0:  
        id, train_ids, labels_dict = write_events_to_file(id, train_ids, labels_dict, train_labels_to_write, train_data_to_write, output_folder)
        
    if len(val_data_to_write) != 0:
        id, val_ids, labels_dict = write_events_to_file(id, val_ids, labels_dict, val_labels_to_write, val_data_to_write, output_folder)
    
    if len(test_data_to_write) != 0:
        id, test_ids, labels_dict = write_events_to_file(id, test_ids, labels_dict, test_labels_to_write, test_data_to_write, output_folder)

    f = open(output_folder+"/labels_dict_"+distribution+".csv", "a")
    writer = csv.writer(f)
    for key, value in labels_dict.items():
        writer.writerow([key] + value)
    f.close()

    f = open(output_folder+"/train_ids_"+distribution+".csv", "a")
    writer = csv.writer(f)
    writer.writerows([[id] for id in train_ids])
    f.close()

    f = open(output_folder+"/test_ids_"+distribution+".csv", "a")
    writer = csv.writer(f)
    writer.writerows([[id] for id in test_ids])
    f.close()

    f = open(output_folder+"/val_ids_"+distribution+".csv", "a")
    writer = csv.writer(f)
    writer.writerows([[id] for id in val_ids])
    f.close()


if __name__ == "__main__":
    """ Extract ttH and ttZ for training of particle assignment. """

    data_folders = ["ttH", "ttZ"]
    productions = ["ttH", "ttZ"]

    for i in range(len(data_folders)):
        data_folder = data_folders[i]
        prod = productions[i]
        file_names = []
        file_names += [join(data_folder, f) for f in listdir(data_folder) if isfile(join(data_folder, f))]

        counter = 1
        for file_name in file_names:
            extract_data(file_name, counter, "data/particle_assignment_training_data", prod, True)
            counter += 1
    
    """ Extract ttW and tt that will then be processed by trained particle
    assignment and the result will be used for mass reco. """
    
    data_folders = ["ttW", "tt"]
    productions = ["ttW", "tt"]

    for i in range(len(data_folders)):
        data_folder = data_folders[i]
        prod = productions[i]
        file_names = []
        file_names += [join(data_folder, f) for f in listdir(data_folder) if isfile(join(data_folder, f))]

        counter = 1
        for file_name in file_names:
            extract_data(file_name, counter, "data/particle_assignment_data_to_be_processed_ttW_tt", prod, False)
            counter += 1
        
    """ Extract ttH and ttZ that will then be processed by trained particle
    assignment and the result will be used for mass reco. """

    data_folders = ["ttH", "ttZ"]
    productions = ["ttH", "ttZ"]

    for i in range(len(data_folders)):
        data_folder = data_folders[i]
        prod = productions[i]
        file_names = []
        file_names += [join(data_folder, f) for f in listdir(data_folder) if isfile(join(data_folder, f))]
    
        counter = 1
        for file_name in file_names:
            extract_data(file_name, counter, "data/particle_assignment_data_to_be_processed_ttH_ttZ", prod, False)
            counter += 1
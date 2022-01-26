/* 
    How to run:
        1) Has to be located in the same folder as folder DiTauMassTools which contains MMC implementation code.
        2) In terminal move to DiTauMassTools/Root and run root. 
        3) In ROOT terminal run following:
            .include /some_path/DiTauMassTools/Root 
            .include /some_path/DiTauMassTools/DiTauMassTools
            .x MissingMassCalculator.cxx
            .x /some_path/mmc.cxx
        4) Calculated mass values are saved to a file.
*/

#include "DiTauMassTools/DiTauMassTools/MissingMassCalculator.h"

#include <iostream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "TH1.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"

// source: https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/
std::vector<std::vector<float_t>> read_csv(std::string filename){

    std::vector<std::vector<float_t>> result;

    // Create an input filestream
    std::cout << filename << std::endl;
    std::ifstream myFile(filename);
    std::cout << filename << std::endl;
    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");
    std::cout << filename << std::endl;
    // Helper vars
    std::string line, colname;
    float_t val;

    int rowIdx = 0;

    // Read data, line by line
    while(std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);
        
        // Extract each integer
        result.push_back(std::vector<float_t> {});

        while(ss >> val){
            
            result.at(rowIdx).push_back(val);
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
            
        }
        rowIdx++;
    }

    // Close file
    myFile.close();

    return result;
}

int numJets25(std::vector<TLorentzVector> jets) {
	int numJets = 0;
	
	TLorentzVector jet;
	for(int i=0; i<jets.size(); ++i) {
		jet = jets[i];
		
		if (jet.E() < 1.0e+07 && jet.Px() < 1.0e+07 && jet.Py() < 1.0e+07 && jet.Pz() < 1.0e+07) {
			if (jet.Et() > 25) {
				numJets++;
			}
		}
			
		
	}
	return numJets;
}

int findIndex(std::vector<int> vec, const int val) {
    int i = 0;
    for(const auto& value: vec) {
        if (value == val) {
            return i;
        }
        i++;
    }
    return -1;
}



// Main function. 
int mmc() {
    
    std::string hello("Hello World!");
    std::cout << hello << std::endl;
	

    MissingMassCalculator MMC;

    std::vector<std::string> output_files = {   "data/mmc/MMC_output_test_narrow_ttH.csv",
                                                "data/mmc/MMC_output_test_narrow_ttZ.csv",
                                                "data/mmc/MMC_output_test_wide_ttH.csv",
                                                "data/mmc/MMC_output_test_wide_ttZ.csv",
                                                "data/mmc/MMC_output_test_wide_ttW.csv",
                                                "data/mmc/MMC_output_test_wide_tt.csv",};

    std::vector<std::string> input_files = {"data/mass_reco/mass_reco_input_narrow_selection_test_ttH.csv",
                                            "data/mass_reco/mass_reco_input_narrow_selection_test_ttZ.csv",
                                            "data/mass_reco/mass_reco_input_wide_selection_test_ttH.csv",
                                            "data/mass_reco/mass_reco_input_wide_selection_test_ttZ.csv",
                                            "data/mass_reco/mass_reco_input_wide_selection_test_ttW.csv",
                                            "data/mass_reco/mass_reco_input_wide_selection_test_tt.csv",};

    for (int i=0; i<6; i++) {
        std::string input_file_name = input_files.at(i);
        std::string output_file_name = output_files.at(i);

        std::ofstream output_file;
        output_file.open(output_file_name);

        std::vector<std::vector<float_t>> data = read_csv(input_file_name);

        for (auto event:data) {

            /* Get the necessary variables. */

            Float_t met_x = event.at(65)/1000;
            Float_t met_y = event.at(66)/1000;

            Float_t tau_px = event.at(28)/1000;
            Float_t tau_py = event.at(29)/1000;
            Float_t tau_pz = event.at(30)/1000;
            Float_t tau_e = event.at(31)/1000;

            Float_t lep_px = event.at(12)/1000;
            Float_t lep_py = event.at(13)/1000;
            Float_t lep_pz = event.at(14)/1000;
            Float_t lep_e = event.at(15)/1000;

            Char_t tau_type = event.at(70);
            Float_t HT = event.at(69)/1000;
            Float_t JetsET = event.at(67)/1000;
            Int_t numJets25 = event.at(68);

            TVector2 met_vec(met_x, met_y);

            TLorentzVector tau, lep;
            tau.SetPxPyPzE(tau_px, tau_py, tau_pz, tau_e);
            lep.SetPxPyPzE(lep_px, lep_py, lep_pz, lep_e);

            Char_t lep_type = 0;
            
            if (lep.Mag() < 6.67/1000) {		// threshold between electron and muon mass
                lep_type = 0;
            } else {
                lep_type = 1;
            }

            /* Following implementation follows the README file located with the MMC implementation. */

            MMC.SetUseEfficiencyRecovery(1);      // to recover efficiency loss 
            MMC.SetCalibrationSet(MMCCalibrationSet::MMC2015);
            
            MMC.SetMetVec(met_vec); // passing MET vector (in form of TVector2 object) 

            MMC.SetVisTauVec(0,tau); // passing TLorentzVec for visible tau-0 (first visible tau in the event) 
            MMC.SetVisTauVec(1,lep); // passing TLorentzVec for visible tau-1 (second visible tau in the event)

            MMC.SetVisTauType(0,tau_type); // passing decay type of tau-0 
            MMC.SetVisTauType(1,lep_type); // passing decay type of tau-1

            Double_t SumEt = HT - JetsET -  tau.Et() - lep.Et();
            MMC.SetSumEt(SumEt);                  // passing event sumEt (RefFinal_BDT_medium, see comments above)
            MMC.SetNjet25(numJets25);                // For Lep-Had analysis only

            int misMassTest=MMC.RunMissingMassCalculator(); // run MMC
            int output_fitstatus=MMC.GetFitStatus(); // MMC output: 1=found solution; 0= no slution

            Double_t MMC_mass;
            TLorentzVector MMC_rec4vec;
            TLorentzVector nu0;
            TLorentzVector nu1;

            MMC_mass = MMC.GetFittedMass(MMCFitMethod::MLM); // returns mass according to default method hh2013 now using meth 1 instead of 0
            MMC_rec4vec = MMC.GetResonanceVec(MMCFitMethod::MAXW); // Optional: returns resonance 4-vec, if you need it
            nu0=MMC.GetNeutrino4vec(MMCFitMethod::MAXW,0); // Optional: returns 4-vec for neutrino-1
            nu1=MMC.GetNeutrino4vec(MMCFitMethod::MAXW,1);
            std::cout << MMC_mass << std::endl;
            output_file << MMC_mass;
            output_file << "\n";

        }

        output_file.close();
    }    

    std::string bye("Bye World!");
    std::cout << bye << std::endl;
    return 0;
}
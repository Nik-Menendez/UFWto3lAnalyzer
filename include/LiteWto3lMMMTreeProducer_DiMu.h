#ifndef LiteWto3lMMMTreeProducer_DiMu_h
#define LiteWto3lMMMTreeProducer_DiMu_h

#include "LiteWto3lMMMTreeProducer_DiMu_Linkdef.h"
#include "Analyzer.h"
#include "HZZTree_fakerate.h"
#include "deltaR.h"
#include <cmath>

using namespace std;  

class LiteWto3lMMMTreeProducer_DiMu : public Analyzer 
{
    public:
        // static const double Zmass = 91.1876;
        const double Zmass = 91.1876;

        LiteWto3lMMMTreeProducer_DiMu();
        ~LiteWto3lMMMTreeProducer_DiMu();
        LiteWto3lMMMTreeProducer_DiMu(
                double isoCutEl_in,
                double isoCutMu_in,
                TString outputDir_in,
                TString outFileName_in,
                bool do_wrong_fc_in=false
                );
         LiteWto3lMMMTreeProducer_DiMu(
                TString outputDir_in,
                TString outFileName_in
                );

        int process();
        void setup();
        void end();
        void initTree();
        bool passSelection();
        void setDebugMode(bool debug_in);
        void sortedArray(double x, double y, double z, double sortarray[3]);

        double isoCutEl=0.35;
        double isoCutMu=0.35;
        double leadingPtCut=20.0; 
        double subleadingPtCut=10.0; 
        double lowestPtCut=5.0;

        TString treeName = "Ana/passedEvents";
        TString outTreeName = "passedEvents";
        bool do_wrong_fc = false;

	double total_events=0;
        double cut1=0;
	double cut2=0;
	double cut3=0;
	double cut4=0;
	double cut5=0;

        //=================== DEBUG ============================
        bool debug=true;
	double what=0;
        double what1=0;
        double what2=0;
	double what3=0;
	double what4=0;
	double huh=0;
	double L3isLoose=0;
        //=================== DEBUG ============================

        TString outFileName;
        TFile* outFile=0;
        TTree* outTree=0;
};

LiteWto3lMMMTreeProducer_DiMu::LiteWto3lMMMTreeProducer_DiMu(){

}

LiteWto3lMMMTreeProducer_DiMu::~LiteWto3lMMMTreeProducer_DiMu(){

}

bool LiteWto3lMMMTreeProducer_DiMu::passSelection(){
    return true;
}


void LiteWto3lMMMTreeProducer_DiMu::initTree(){
    setHZZTree(tree);
}

void LiteWto3lMMMTreeProducer_DiMu::setDebugMode(bool debug_in){
    debug = debug_in;
}

void LiteWto3lMMMTreeProducer_DiMu::setup(){
    outFile = TFile::Open(outputDir+outFileName,fOptionWrite);
    outTree = new TTree(outTreeName,outTreeName);

    initNewLiteTree_fakerate(outTree);
}

void LiteWto3lMMMTreeProducer_DiMu::end(){
    outFile->cd();
    outTree->Write(outTreeName,TObject::kOverwrite);
    outFile->Close(); 

    double pass1 = total_events-cut1;
    double pass2 = pass1-cut2;
    double pass3 = pass2-cut3;
    double pass4 = pass3-cut4;
    double pass5 = pass4-cut5;

    if(!debug){
    std::cout << std::endl;
    std::cout << "Efficiencies for each cut" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "Total events before cuts: " << total_events << std::endl;
    if (total_events!=0) std::cout << "Pass 2 lepton cut: " << pass1 << ". Efficiency = " << (pass1/total_events)*100 << "%" << std::endl;
    if (pass1!=0) std::cout << "Pass trigger: " << pass2 << ". Efficiency = " << (pass2/pass1)*100 << "%" << std::endl;
    if (pass2!=0) std::cout << "Pass tight muon cut: " << pass3 << ". Efficiency = " << (pass3/pass2)*100 << "%" << std::endl;
    if (pass3!=0) std::cout << "Pass muons different charge cut: " << pass4 << ". Efficiency = " << (pass4/pass3)*100 << "%" << std::endl;
    if (pass4!=0) std::cout << "Pass muon pt cut: " << pass5 << ". Efficiency = " << (pass5/pass4)*100 << "%" << std::endl;
    std::cout << "=====================================================" << std::endl;
    if (total_events!=0) std::cout << "Total Efficiency = " << (pass5/total_events)*100 << "%" << std::endl; 
    std::cout << "Percent of events where lep3 has Iso > 0.35: " << (L3isLoose/pass5)*100 << "%" << std::endl;
    } else {
    //=================== DEBUG ============================
    std::cout << std::endl;
    std::cout << (what/total_events)*100 << "% of events had different sized id and pt???" << std::endl;
    std::cout << (what1/what)*100 << "% of the time id is bigger" << std::endl;
    std::cout << (what2/what)*100 << "% of the time pt is bigger" << std::endl;
    std::cout << (what3/what)*100 << "% of the time id and pt differ by 3" << std::endl;
    std::cout << (what4/what)*100 << "% of the time id and pt differ by more than 3" << std::endl;
    std::cout << (huh/total_events)*100 << "% of events have no leptons" << std::endl;
    //=================== DEBUG ============================
    }
}

LiteWto3lMMMTreeProducer_DiMu::LiteWto3lMMMTreeProducer_DiMu(
                double isoCutEl_in,
                double isoCutMu_in,
                TString outputDir_in,
                TString outFileName_in,
                bool do_wrong_fc_in
                ){
    isoCutEl    = isoCutEl_in;
    isoCutMu    = isoCutMu_in;
    outputDir   = outputDir_in;
    outFileName = outFileName_in;
    do_wrong_fc = do_wrong_fc_in;
}

LiteWto3lMMMTreeProducer_DiMu::LiteWto3lMMMTreeProducer_DiMu(
                TString outputDir_in,
                TString outFileName_in
                ){
    outputDir   = outputDir_in;
    outFileName = outFileName_in;
}

void LiteWto3lMMMTreeProducer_DiMu::sortedArray(double x, double y, double z, double sortarray[3])
{
    double max_v = max(x,max(y,z));
    double min_v = min(x,min(y,z));
    double mid_v = -1;
    if (x != max_v && x != min_v) mid_v = x;
    if (y != max_v && y != min_v) mid_v = y;
    if (z != max_v && z != min_v) mid_v = z;

    sortarray[0] = max_v;
    sortarray[1] = mid_v;
    sortarray[2] = min_v;
    //return sortarray;
}

int LiteWto3lMMMTreeProducer_DiMu::process(){

    total_events++;

    int Nlep = (*lep_id).size();
    int Nlep_diff = (*lep_id).size() - (*lep_pt).size();
    int nTightLep = 0;
    int nLooseLep = 0;
    vector<int> tightIsoLepIndex;
    vector<int> looseIsoLepIndex;

    nLeptons = Nlep;

    //=================== DEBUG ============================

    if(debug){
    if(((*lep_id).size() != (*lep_pt).size())){
        what++;
        if((*lep_id).size() > (*lep_pt).size()){what1++;}
        else{what2++;}
        if(((*lep_id).size() - (*lep_pt).size()) == 3){what3++;}
        else if(((*lep_id).size() - (*lep_pt).size()) > 3){what4++;}
        //std::cout << "lep_id.size() = " << (*lep_id).size() << ", lep_pt.size() = " << (*lep_pt).size() << ", difference = " << Nlep_diff << std::endl; 
    }
    nLeptons_id = (*lep_id).size();
    nLeptons_pt = (*lep_pt).size();
    nLeptons_diff = Nlep_diff;
    if((*lep_id).size()==0){huh++;}

    if((*lep_id).size()>=3){
        nLeptons_diff_pass = (*lep_pt).size();
    } else {
        nLeptons_diff_pass = -1;
    }
    }

    //=================== DEBUG ============================


    if(!debug){
    if((*lep_id).size()<2 || (*lep_pt).size()<2){cut1++;return -1;}
    //if(Nlep<3){cut1++;return -1;}

    if(passedTrig == 0){cut2++;return -1;}

    for (int iLep = 0; iLep < Nlep; iLep++) {
        if ((*lep_tightId)[iLep] == 1 && (*lep_RelIso)[iLep] < 0.35 && abs((*lep_id)[iLep]) == 13/* && nTightLep < 3*/){
            nTightLep++;
            tightIsoLepIndex.push_back(iLep);
        }
        //else if ((*lep_tightId)[iLep] == 1 && abs((*lep_id)[iLep]) == 13){
        //    nLooseLep++;
        //    looseIsoLepIndex.push_back(iLep);
        //}
    }

    if (!(nTightLep >= 2/* && nLooseLep >= 1*/)){cut3++;return -1;}

    int index1 = tightIsoLepIndex[0];
    int index2 = tightIsoLepIndex[1];
	//int index3 = 0;
    //if (nTightLep >= 3){index3 = tightIsoLepIndex[2];}

    //if ( (*lep_id)[index1] + (*lep_id)[index2] != 0 ) return -1;

	//if (nTightLep >= 3) {
    //	if ((*lep_id)[index1] > 0 && (*lep_id)[index2] > 0 && (*lep_id)[index3] > 0){cut4++;return -1;}
    //	if ((*lep_id)[index1] < 0 && (*lep_id)[index2] < 0 && (*lep_id)[index3] < 0){cut4++;return -1;}
	//} else {
		if ((*lep_id)[index1] > 0 && (*lep_id)[index2] > 0){cut4++;return -1;}
		if ((*lep_id)[index1] < 0 && (*lep_id)[index2] < 0){cut4++;return -1;}
	//}

	//if (nTightLep >= 3) {
	//	if ((*lep_id)[index1] > 0 && (*lep_id)[index2] > 0){index2 = index3;}
	//	if ((*lep_id)[index1] < 0 && (*lep_id)[index2] < 0){index2 = index3;}
	//}

    //if ((*lep_Sip)[index1] > 3 || (*lep_Sip)[index2] > 3 || (*lep_Sip)[index3] > 3) return -1;

    //double pTs[3]; sortedArray((*lep_pt)[index1], (*lep_pt)[index2], (*lep_pt)[index3], pTs);
    //if (pTs[0] < leadingPtCut || pTs[1] < subleadingPtCut || pTs[2] < lowestPtCut){cut5++;return -1;}
	if ((*lep_pt)[index1] < leadingPtCut || (*lep_pt)[index2] < subleadingPtCut){cut5++;return -1;}

    //if((*lep_RelIso)[index3]>0.35){L3isLoose++;}


    //double isos[3]; sortedArray((*lep_RelIso)[index1], (*lep_RelIso)[index2], (*lep_RelIso)[index3], isos);
    //double sips[3]; sortedArray((*lep_Sip)[index1], (*lep_Sip)[index2], (*lep_Sip)[index3], sips);

    int tmp = 0;
    TLorentzVector Lep1,Lep2,Lep3;
    //Lep3.SetPtEtaPhiM((*lep_pt)[index3],(*lep_eta)[index3],(*lep_phi)[index3],(*lep_mass)[index3]);
    if ((*lep_pt)[index1] >= (*lep_pt)[index2]){
        Lep1.SetPtEtaPhiM((*lep_pt)[index1],(*lep_eta)[index1],(*lep_phi)[index1],(*lep_mass)[index1]); 
        Lep2.SetPtEtaPhiM((*lep_pt)[index2],(*lep_eta)[index2],(*lep_phi)[index2],(*lep_mass)[index2]);
    }
    if ((*lep_pt)[index1] < (*lep_pt)[index2]){
        tmp = index1;
        index1 = index2;
        index2 = tmp;
        Lep1.SetPtEtaPhiM((*lep_pt)[index1],(*lep_eta)[index1],(*lep_phi)[index1],(*lep_mass)[index1]); 
        Lep2.SetPtEtaPhiM((*lep_pt)[index2],(*lep_eta)[index2],(*lep_phi)[index2],(*lep_mass)[index2]);
    }

    TLorentzVector Leps = Lep1+Lep2+Lep3;
    double LepsPhi = Leps.Phi(); double LepsPt = Leps.Pt();
    //double mT = sqrt(2*(met)*LepsPt*(1-cos(deltaPhi(LepsPhi, met_phi) ) ) );

    vector<double> dRs; double tmpDr = -1;
    vector<TLorentzVector> mlls;
    if ((*lep_id)[index1] != (*lep_id)[index2]) {
        mlls.push_back(Lep1 + Lep2);
        dRs.push_back(deltaR(Lep1.Eta(),Lep1.Phi(),Lep2.Eta(),Lep2.Phi()));

    } else {tmpDr = deltaR(Lep1.Eta(),Lep1.Phi(),Lep2.Eta(),Lep2.Phi());}

    //if ((*lep_id)[index1] != (*lep_id)[index3]) {
    //    mlls.push_back(Lep1 + Lep3);
    //    dRs.push_back(deltaR(Lep1.Eta(),Lep1.Phi(),Lep3.Eta(),Lep3.Phi()));
  
    //} else {tmpDr = deltaR(Lep1.Eta(),Lep1.Phi(),Lep3.Eta(),Lep3.Phi());}

    //if ((*lep_id)[index2] != (*lep_id)[index3]) {
    //    mlls.push_back(Lep2 + Lep3);
    //    dRs.push_back(deltaR(Lep2.Eta(),Lep2.Phi(),Lep3.Eta(),Lep3.Phi()));

    //} else {tmpDr = deltaR(Lep2.Eta(),Lep2.Phi(),Lep3.Eta(),Lep3.Phi());}

    dRs.push_back(tmpDr);

    TLorentzVector Z = Lep1+Lep2;
    //if (!(Z.M() > 86 && Z.M() < 96)) return kTRUE;

    idL1 = (*lep_id)[index1]; pTL1 = Lep1.Pt(); etaL1 = Lep1.Eta();
    idL2 = (*lep_id)[index2]; pTL2 = Lep2.Pt(); etaL2 = Lep2.Eta();
    //idL3 = (*lep_id)[index3]; pTL3 = Lep3.Pt(); etaL3 = Lep3.Eta();       
    phiL1 = Lep1.Phi();//deltaphiL13 = deltaPhi (Lep1.Phi(), Lep3.Phi()); 
    phiL2 = Lep2.Phi();//deltaphiL14 = deltaPhi (Lep1.Phi(), Lep4.Phi());
    phiL3 = Lep3.Phi();//deltaphiL23 = deltaPhi (Lep2.Phi(), Lep3.Phi());
    vector<TLorentzVector> P4s; vector<int> tmpIDs;             
    P4s.push_back(Lep1); P4s.push_back(Lep2);
    P4s.push_back(Lep3); 
    tmpIDs.push_back(idL1); tmpIDs.push_back(idL2);
    tmpIDs.push_back(idL3); 
    IsoL1 = (*lep_RelIso)[index1]; IsoL2 = (*lep_RelIso)[index2];
    //IsoL3 = (*lep_RelIso)[index3];
    massL1 = (*lep_mass)[index1]; massL2 = (*lep_mass)[index2];
    //massL3 = (*lep_mass)[index3];
    MomIdL1 = (*lep_matchedR03_MomId)[index1];  MomIdL2 = (*lep_matchedR03_MomId)[index2];  //MomIdL3 = (*lep_matchedR03_MomId)[index3];
    PDG_IdL1 = (*lep_matchedR03_PdgId)[index1];  PDG_IdL2 = (*lep_matchedR03_PdgId)[index2];  //PDG_IdL3 = (*lep_matchedR03_PdgId)[index3];
    MomMomIdL1 = (*lep_matchedR03_MomMomId)[index1];  MomMomIdL2 = (*lep_matchedR03_MomMomId)[index2];  //MomMomIdL3 = (*lep_matchedR03_MomMomId)[index3];
    dR12 = deltaR(Lep1.Eta(),Lep1.Phi(),Lep2.Eta(),Lep2.Phi()); dR13 = deltaR(Lep1.Eta(),Lep1.Phi(),Lep3.Eta(),Lep3.Phi()); dR23 = deltaR(Lep2.Eta(),Lep2.Phi(),Lep3.Eta(),Lep3.Phi());

    TLorentzVector threeleps = Lep1+Lep2+Lep3;
    m3l = threeleps.M();

    TLorentzVector Met;
    Met.SetPtEtaPhiM(met,0,met_phi,0);
    mt = (threeleps + Met).M();

    massZ1 = Z.M();
    pT3l = Leps.Pt();

    }
           
    //cout<<"fill tree"<<endl;
    //cout<<endl;
    outTree->Fill();   
return 1;
}

#endif

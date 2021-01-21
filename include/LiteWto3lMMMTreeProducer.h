#ifndef LiteWto3lMMMTreeProducer_h
#define LiteWto3lMMMTreeProducer_h

#include "LiteWto3lMMMTreeProducer_Linkdef.h"
#include "Analyzer.h"
#include "HZZTree_fakerate.h"
#include "deltaR.h"
#include <cmath>

using namespace std;  

class LiteWto3lMMMTreeProducer : public Analyzer 
{
    public:
        // static const double Zmass = 91.1876;
        const double Zmass = 91.1876;

        LiteWto3lMMMTreeProducer();
        ~LiteWto3lMMMTreeProducer();
        LiteWto3lMMMTreeProducer(
                double isoCutEl_in,
                double isoCutMu_in,
                TString outputDir_in,
                TString outFileName_in,
                bool do_wrong_fc_in=false
                );
         LiteWto3lMMMTreeProducer(
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

LiteWto3lMMMTreeProducer::LiteWto3lMMMTreeProducer(){

}

LiteWto3lMMMTreeProducer::~LiteWto3lMMMTreeProducer(){

}

bool LiteWto3lMMMTreeProducer::passSelection(){
    return true;
}


void LiteWto3lMMMTreeProducer::initTree(){
    setHZZTree(tree);
}

void LiteWto3lMMMTreeProducer::setDebugMode(bool debug_in){
    debug = debug_in;
}

void LiteWto3lMMMTreeProducer::setup(){
    outFile = TFile::Open(outputDir+outFileName,fOptionWrite);
    outTree = new TTree(outTreeName,outTreeName);

    initNewLiteTree_fakerate(outTree);
}

void LiteWto3lMMMTreeProducer::end(){
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
    if (total_events!=0) std::cout << "Pass 3 lepton cut: " << pass1 << ". Efficiency = " << (pass1/total_events)*100 << "%" << std::endl;
    if (pass1!=0) std::cout << "Pass trigger: " << pass2 << ". Efficiency = " << (pass2/pass1)*100 << "%" << std::endl;
    if (pass2!=0) std::cout << "Pass tight lepton cut: " << pass3 << ". Efficiency = " << (pass3/pass2)*100 << "%" << std::endl;
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

LiteWto3lMMMTreeProducer::LiteWto3lMMMTreeProducer(
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

LiteWto3lMMMTreeProducer::LiteWto3lMMMTreeProducer(
                TString outputDir_in,
                TString outFileName_in
                ){
    outputDir   = outputDir_in;
    outFileName = outFileName_in;
}

void LiteWto3lMMMTreeProducer::sortedArray(double x, double y, double z, double sortarray[3])
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

int LiteWto3lMMMTreeProducer::process(){

    total_events++;

    unsigned int Nlep = (*lep_id).size();
    int Nlep_diff = (*lep_id).size() - (*lep_pt).size();
    int nTightLep = 0;
    int nLooseLep = 0;
    vector<int> tightIsoLepIndex;
    vector<int> looseIsoLepIndex;

    nLeptons = Nlep;

	int index1=0; int index2=0; int index3=0;
	bool foundZ1LCandidate=false;

	//*******************************************

	const double Zmass = 91.1876;
	double mZ1Low = 0.0;//40.0;
	double mZ1High = 200.0;//120.0;
	bool onlymu = false;

	if ( Nlep<3 ) {cut1++; return -1;}

	int n_Zs=0;
	vector<int> Z_Z1L_lepindex1;
	vector<int> Z_Z1L_lepindex2;

	for(unsigned int i=0; i<Nlep; i++){
        for(unsigned int j=i+1; j<Nlep; j++){

            // same flavor opposite charge
            if(((*lep_id)[i]+(*lep_id)[j])!=0) continue;
			if(onlymu) {
				if(abs((*lep_id)[i])!=13) continue;
			}

            TLorentzVector li, lj;
            li.SetPtEtaPhiM((*lep_pt)[i],(*lep_eta)[i],(*lep_phi)[i],(*lep_mass)[i]);
            lj.SetPtEtaPhiM((*lep_pt)[j],(*lep_eta)[j],(*lep_phi)[j],(*lep_mass)[j]);

            TLorentzVector Z;
            Z = li+lj;

            if (Z.M()>0.0) {
                n_Zs++;
                Z_Z1L_lepindex1.push_back(i);
                Z_Z1L_lepindex2.push_back(j);
            }

        } // lep i
    } // lep j

	bool properLep_ID = false; int Nmm = 0; int Nmp = 0; int Nem = 0; int Nep = 0;
    for(unsigned int i =0; i<Nlep; i++) {
		if(abs((*lep_id)[i])!=13) continue;
        if((*lep_id)[i]<0) Nmm = Nmm+1;
        if((*lep_id)[i]>0) Nmp = Nmp+1;
    }
	for(unsigned int i =0; i<Nlep; i++) {
        if(abs((*lep_id)[i])!=11) continue;
        if((*lep_id)[i]<0) Nem = Nem+1;
        if((*lep_id)[i]>0) Nep = Nep+1;
    }

	if(Nmm>=1 && Nmp>=1) properLep_ID = true; //2mu + x
	if(!onlymu) {
		if(Nem>=1 && Nep>=1) properLep_ID = true; //2e + x
	}

	if(!properLep_ID) {cut4++; return -1;}

	// Consider all Z candidates
    double minZ1DeltaM=9999.9;
    for (int i=0; i<n_Zs; i++) {

        int i1 = Z_Z1L_lepindex1[i]; int i2 = Z_Z1L_lepindex2[i];
		int j1 = 3 - i1 - i2;
		if(onlymu){
			if(abs((*lep_id)[j1])!=13) continue;
		}
        //int j1 = 999; //3 - i1 - i2; // index of the third lepton (check if this works)
		//for (unsigned int j=0;j<Nlep;j++) {
		//	if(Nlep==i1 || Nlep==i2) continue;
		//	if(abs((*lep_id)[j])==13) {
		//		j1 = Nlep;
		//		break;
		//	}
		//}
		//if(j1==999) continue;

        TLorentzVector lep_i1, lep_i2, lep_j1;
        lep_i1.SetPtEtaPhiM((*lep_pt)[i1],(*lep_eta)[i1],(*lep_phi)[i1],(*lep_mass)[i1]);
        lep_i2.SetPtEtaPhiM((*lep_pt)[i2],(*lep_eta)[i2],(*lep_phi)[i2],(*lep_mass)[i2]);
        lep_j1.SetPtEtaPhiM((*lep_pt)[j1],(*lep_eta)[j1],(*lep_phi)[j1],(*lep_mass)[j1]);

        TLorentzVector Zi;
        Zi = lep_i1+lep_i2;
        //Zi.SetPtEtaPhiM(Z_Z1L_pt[i],Z_Z1L_eta[i],Z_Z1L_phi[i],Z_Z1L_mass[i]);

        TLorentzVector Z1 = Zi;
        double Z1DeltaM = abs(Zi.M()-Zmass);
        int Z1_lepindex[2] = {0,0};
        if (lep_i1.Pt()>lep_i2.Pt()) { Z1_lepindex[0] = i1;  Z1_lepindex[1] = i2; }
        else { Z1_lepindex[0] = i2;  Z1_lepindex[1] = i1; }

        // Check Leading and Subleading pt Cut
        vector<double> allPt;
        allPt.push_back(lep_i1.Pt()); allPt.push_back(lep_i2.Pt());
        std::sort(allPt.begin(), allPt.end());
        if (debug) cout<<" leading pt: "<<allPt[1]<<" cut: "<<leadingPtCut<<" subleadingPt: "<<allPt[0]<<" cut: "<<subleadingPtCut<<endl;
        if (allPt[1]<leadingPtCut || allPt[0]<subleadingPtCut ) continue;

        // Check dR(li,lj)>0.02 for any i,j
        vector<double> alldR;
        alldR.push_back(deltaR(lep_i1.Eta(),lep_i1.Phi(),lep_i2.Eta(),lep_i2.Phi()));
        alldR.push_back(deltaR(lep_i1.Eta(),lep_i1.Phi(),lep_j1.Eta(),lep_j1.Phi()));
        alldR.push_back(deltaR(lep_i2.Eta(),lep_i2.Phi(),lep_j1.Eta(),lep_j1.Phi()));
        if (debug) cout<<" minDr: "<<*min_element(alldR.begin(),alldR.end())<<endl;
        if (*min_element(alldR.begin(),alldR.end())<0.02) continue;

        // Check M(l+,l-)>4.0 GeV for any OS pair
        // Do not include FSR photons
        vector<double> allM;
        TLorentzVector i1i2;
        i1i2 = (lep_i1)+(lep_i2); allM.push_back(i1i2.M());
        if ((*lep_id)[i1]*(*lep_id)[j1]<0) {
            TLorentzVector i1j1;
            i1j1 = (lep_i1)+(lep_j1); allM.push_back(i1j1.M());
        } else {
            TLorentzVector i2j1;
            i2j1 = (lep_i2)+(lep_j1); allM.push_back(i2j1.M());
        }
        if (debug) cout<<" min m(l+l-): "<<*min_element(allM.begin(),allM.end())<<endl;
        if (*min_element(allM.begin(),allM.end())<4.0) {continue;}

        // Check isolation cut (without FSR ) for Z1 leptons
        if ((*lep_RelIso)[Z1_lepindex[0]]>((abs((*lep_id)[Z1_lepindex[0]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
        if ((*lep_RelIso)[Z1_lepindex[1]]>((abs((*lep_id)[Z1_lepindex[1]])==11) ? isoCutEl : isoCutMu)) continue; // checking iso with FSR removed
        // Check tight ID cut for Z1 leptons
        if (!((*lep_tightId)[Z1_lepindex[0]])) continue; // checking tight lepton ID
        if (!((*lep_tightId)[Z1_lepindex[1]])) continue; // checking tight lepton ID

        if ( (Z1.M() < mZ1Low) || (Z1.M() > mZ1High) ) continue;

        if (debug) cout<<"good Z1L candidate, Z1DeltaM: "<<Z1DeltaM<<" minZ1DeltaM: "<<minZ1DeltaM<<endl;

        // Check if this candidate has the best Z1 and highest scalar sum of Z2 lepton pt

        if ( Z1DeltaM<=minZ1DeltaM ) {

            minZ1DeltaM = Z1DeltaM;

            TLorentzVector Z1L;
            Z1L = Z1+lep_j1;

            massZ1_Z1L = Z1.M();
            mass3l = Z1L.M();

            index1 = Z1_lepindex[0];
            index2 = Z1_lepindex[1];
            index3 = j1;

            if (debug) cout<<" new best Z1L candidate: massZ1: "<<massZ1<<" (mass3l: "<<mass3l<<")"<<endl;
            foundZ1LCandidate=true;

        }
    }

	//*******************************************

	if(!foundZ1LCandidate) {cut5++; return -1;}

    int tmp = 0;
    TLorentzVector Lep1,Lep2,Lep3;
    Lep3.SetPtEtaPhiM((*lep_pt)[index3],(*lep_eta)[index3],(*lep_phi)[index3],(*lep_mass)[index3]);
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

    if ((*lep_id)[index1] != (*lep_id)[index3]) {
        mlls.push_back(Lep1 + Lep3);
        dRs.push_back(deltaR(Lep1.Eta(),Lep1.Phi(),Lep3.Eta(),Lep3.Phi()));
  
    } else {tmpDr = deltaR(Lep1.Eta(),Lep1.Phi(),Lep3.Eta(),Lep3.Phi());}

    if ((*lep_id)[index2] != (*lep_id)[index3]) {
        mlls.push_back(Lep2 + Lep3);
        dRs.push_back(deltaR(Lep2.Eta(),Lep2.Phi(),Lep3.Eta(),Lep3.Phi()));

    } else {tmpDr = deltaR(Lep2.Eta(),Lep2.Phi(),Lep3.Eta(),Lep3.Phi());}

    dRs.push_back(tmpDr);

    TLorentzVector Z = Lep1+Lep2;
    //if (!(Z.M() > 86 && Z.M() < 96)) return kTRUE;

    idL1 = (*lep_id)[index1]; pTL1 = Lep1.Pt(); etaL1 = Lep1.Eta();
    idL2 = (*lep_id)[index2]; pTL2 = Lep2.Pt(); etaL2 = Lep2.Eta();
    idL3 = (*lep_id)[index3]; pTL3 = Lep3.Pt(); etaL3 = Lep3.Eta();       
    phiL1 = Lep1.Phi();//deltaphiL13 = deltaPhi (Lep1.Phi(), Lep3.Phi()); 
    phiL2 = Lep2.Phi();//deltaphiL14 = deltaPhi (Lep1.Phi(), Lep4.Phi());
    phiL3 = Lep3.Phi();//deltaphiL23 = deltaPhi (Lep2.Phi(), Lep3.Phi());
    vector<TLorentzVector> P4s; vector<int> tmpIDs;             
    P4s.push_back(Lep1); P4s.push_back(Lep2);
    P4s.push_back(Lep3); 
    tmpIDs.push_back(idL1); tmpIDs.push_back(idL2);
    tmpIDs.push_back(idL3); 
    IsoL1 = (*lep_RelIso)[index1]; IsoL2 = (*lep_RelIso)[index2];
    IsoL3 = (*lep_RelIso)[index3];
	tightIdL1 = (*lep_tightId)[index1]; tightIdL2 = (*lep_tightId)[index2]; tightIdL3 = (*lep_tightId)[index3];
    massL1 = (*lep_mass)[index1]; massL2 = (*lep_mass)[index2];
    massL3 = (*lep_mass)[index3];
    MomIdL1 = (*lep_matchedR03_MomId)[index1];  MomIdL2 = (*lep_matchedR03_MomId)[index2];  MomIdL3 = (*lep_matchedR03_MomId)[index3];
    PDG_IdL1 = (*lep_matchedR03_PdgId)[index1];  PDG_IdL2 = (*lep_matchedR03_PdgId)[index2];  PDG_IdL3 = (*lep_matchedR03_PdgId)[index3];
    MomMomIdL1 = (*lep_matchedR03_MomMomId)[index1];  MomMomIdL2 = (*lep_matchedR03_MomMomId)[index2];  MomMomIdL3 = (*lep_matchedR03_MomMomId)[index3];
    dR12 = deltaR(Lep1.Eta(),Lep1.Phi(),Lep2.Eta(),Lep2.Phi()); dR13 = deltaR(Lep1.Eta(),Lep1.Phi(),Lep3.Eta(),Lep3.Phi()); dR23 = deltaR(Lep2.Eta(),Lep2.Phi(),Lep3.Eta(),Lep3.Phi());

    TLorentzVector threeleps = Lep1+Lep2+Lep3;
    m3l = threeleps.M();

    TLorentzVector Met;
    Met.SetPtEtaPhiM(met,0,met_phi,0);
    mt = (threeleps + Met).M();

    massZ1 = Z.M();
    pT3l = Leps.Pt();

    //cout<<"fill tree"<<endl;
    //cout<<endl;
    outTree->Fill();   
return 1;
}

#endif

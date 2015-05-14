#include <iostream>
#include <iomanip>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

void trackExtraction(const char *file) {
    
    TFile *f = new TFile(file);
    TTree *tree = (TTree*) f->Get("T");
    
    RAT::DS::Root *rds = new RAT::DS::Root();
    tree->SetBranchAddress("ds", &rds);
    
    int nEvents = tree->GetEntries();
    
    for (int i = 0; i < nEvents; i++) {

        tree->GetEntry(i);
        RAT::DS::MC *mc = rds->GetMC();
        RAT::DS::MCParticle *prim = mc->GetMCParticle(0);
 
        printf("MC PE     : %8.0f\n", mc->GetNumPE());
        printf("MC PMT Cnt: %8.0f\n", mc->GetMCPMTCount());
        printf("MC Trk Cnt: %8.0f\n", mc->GetMCTrackCount());
        printf("MCPart Cnt: %8.0f\n", mc->GetMCParticleCount());
        printf("Kin E     : %8.2f\n", prim->GetKE());
        printf("position  : %8.2f %8.2f %8.2f\n", prim->GetPosition().X(), prim->GetPosition().Y(), prim->GetPosition().Z());
        printf("momentum  : %8.3f %8.3f %8.3f\n\n", prim->GetMomentum().X(), prim->GetMomentum().Y(), prim->GetMomentum().Z());
        
    }
}

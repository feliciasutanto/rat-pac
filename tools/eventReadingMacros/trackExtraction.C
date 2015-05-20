// Extract the tracks from a file using root command prompts
// example of use. The command
// >  root trackExtraction.C\(\"/Users/marcbergevin/RAT_ROOT/output.root\"\)
// Will apply this routine of the output.root file
// M.B
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
        
        
        int nTracks = mc->GetMCTrackCount();
        
        //        Particle *trackmap = new Particle[nTracks+1];
        for (int j = 0; j < nTracks; j++) {
            RAT::DS::MCTrack *track = mc->GetMCTrack(j);
            int tid = track->GetID();
            int pid = track->GetParentID();
            
            RAT::DS::MCTrackStep *first = track->GetMCTrackStep(0);
            RAT::DS::MCTrackStep *last = track->GetLastMCTrackStep();
            
            if((track->GetParticleName() != "opticalphoton" ) ){
                if(first->ke>1.0 || track->pdgcode>22){//!((last->GetProcess()=="eIoni")||(last->GetProcess()=="hIoni"))){//
                    printf("%7d  %7d %7d %7d %10d %8.3f %10.3f\n",i,j,pid,tid,track->pdgcode, first->ke,first->globalTime);//
                }
            }// if((track->GetParticleName() != "opticalphoton" ) ){
            int nSteps = track->GetMCTrackStepCount();
            TVector3 *steps = new TVector3[nSteps];
            for (int k = 0; k < nSteps; k++) {
                steps[k] = track->GetMCTrackStep(k)->GetEndpoint();
            }//for (int k = 0; k < nSteps; k++)
            
        }//for (int j = 0; j < nTracks; j++) {
    }
}

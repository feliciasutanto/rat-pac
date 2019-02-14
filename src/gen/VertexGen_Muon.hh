#ifndef __RAT_VertexGen_Muon__
#define __RAT_VertexGen_Muon__

#include <RAT/GLG4VertexGen.hh>
#include "RAT/DB.hh"
//#include <RAT/MuonMessenger.hh>

#include <G4Event.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>
#include <vector>
#include <TF1.h>
#include <TF2.h>

namespace RAT {
    
    class VertexGen_Muon : public GLG4VertexGen {
        
    public:
        VertexGen_Muon( const char *arg_dbname = "muon" );
        virtual 		~VertexGen_Muon();
        virtual void 	GeneratePrimaryVertex( G4Event* argEvent,G4ThreeVector& dx,G4double dt);
        virtual void 	SetState( G4String newValues ); //set the state for generator
        virtual 		G4String GetState(); //return current state
        
    private:
        void CalculateMuonVariables();
        void Initialization( bool igflag, double s_hor, double s_ver1, double s_ver2 );
        void Sampling( double *E, double *theta, double *phi, double *dep );
        
        void SetMuonX0( double tmp ) {m_x0=tmp;};
        double GetMuonX0()           {return m_x0;};
        void SetMuonY0( double tmp ) {m_y0=tmp;};
        double GetMuonY0()           {return m_y0;};
        void SetMuonZ0( double tmp ) {m_z0=tmp;};
        double GetMuonZ0()           {return m_z0;};
        void SetMuonE( double tmp )  {m_E=tmp;};
        double GetMuonE()            {return m_E;};
        void SetMuonSign( int tmp )  {m_signid=tmp;};
        int GetMuonSign()            {return m_signid;};
        void SetMuonCx( double tmp ) {m_cx=tmp;};
        double GetMuonCx()           {return m_cx;};
        void SetMuonCy(double tmp )  {m_cy=tmp;};
        double GetMuonCy()           {return m_cy;};
        void SetMuonCz( double tmp ) {m_cz=tmp;};
        double GetMuonCz()           {return m_cz;};
        
        double m_x0,m_y0,m_z0;
        double m_E,m_signid;
        double m_cx,m_cy,m_cz;
        
        //Global variables
        double spmu[121][62][51], fnmu[32401], depth[360][91], fmu[360][91], mu_cdf[32401],
        e1, e2, m_theta1, m_theta2, m_phi1, m_phi2;
        
        G4ParticleDefinition *ionDef;
        G4ParticleDefinition *neutronDef;
        G4ParticleDefinition *muonplusDef;
        G4ParticleDefinition *muonminusDef;
        G4ParticleDefinition *gammaDef;
        G4ParticleDefinition *n;
        
        //MuonMessenger* messenger;
        G4ThreeVector nu_dir; //useless

    };
    
} // namespace RAT

#endif

#ifndef SPAnalysis_h
#define SPAnalysis_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

// ROOT includes
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TF1.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

/**  This is the header associated with the cc file I posted.
 *   There's not much to it except one little caveat:
 *   Any variables you assign xml-registered parameters to
 *   (in the constructor) must be declared in this file
 *   (look to the very bottom under the "protected" tag). 
 * 
 * @author F. Gaede, DESY
 * @version $Id: BasicAnalysis.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class SPAnalysis : public Processor {

    public:

        virtual Processor*  newProcessor() { return new SPAnalysis ; }


        SPAnalysis() ;

        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init() ;

        /** Called for every run.
        */
        virtual void processRunHeader( LCRunHeader* run ) ;

        /** Called for every event - the working horse.
        */
        virtual void processEvent( LCEvent * evt ) ; 


        virtual void check( LCEvent * evt ) ; 


        /** Called after data processing for clean up.
        */
        virtual void end() ;


    protected:

        #define HCAL_MAX_LAYERS 40
        #define MUON_MAX_LAYERS 12
        #define ECAL_MAX_LAYERS 40

        /** Input collection name.
        */
        std::string _colName ;
        string ROOToutputfilename;
        string ParticleParameters; // particle type and energy
	string Detectorname;
        float Multiplier;
        
        // HCALvsECAL by request for SiD Opt Meeting 03/22/17
        
        TH2F* ECalvsHCalEnergy;
        
        // Collection Names
        string HCALBarrelCollectionName;
        string HCALBarrelRecoCollectionName;
        string HCALEndcapRecoCollectionName;
        string ECALBarrelRecoCollectionName;
        string ECALEndcapRecoCollectionName;
        string HCALEndcapCollectionName;
        string ECALBarrelCollectionName;
        string ECALEndcapCollectionName;
        string MCParticleCollectionName;
        
        TFile* rootoutputfile;
        
        TH1F* DetectedMuonSignals;
        TH1F* DetectedMuonSignalsReco;
        
        // Total Energy Plots
            TH1F* HCalTotalEnergy;
            TH1F* HCalBarrelTotalEnergy;
            TH1F* HCalEndcapTotalEnergy;
            
            TH1F* HCalRecoTotalEnergy;            
            TH1F* HCalEndcapRecoTotalEnergy;
            TH1F* HCalBarrelRecoTotalEnergy;            

            TH1F* ECalTotalEnergy;
            TH1F* ECalBarrelTotalEnergy;
            TH1F* ECalEndcapTotalEnergy;            
            
            TH1F* ECalRecoTotalEnergy;
            TH1F* ECalBarrelRecoTotalEnergy;
            TH1F* ECalEndcapRecoTotalEnergy; 
            
            TH1F* TotalEnergy;
            TH1F* RecoTotalEnergy;
        
        //Energy Distribution Plots
            TH1F* RecoTotalEnergyDistribution;
            TH1F* TotalEnergyDistribution; //TODO
            
            TH1F* HCalEnergyDistribution; //TODO
            TH1F* HCalBarrelEnergyDistribution;
            TH1F* HCalEndcapEnergyDistribution;

            TH1F* HCalRecoEnergyDistribution; //TODO
            TH1F* HCalBarrelRecoEnergyDistribution;            
            TH1F* HCalEndcapRecoEnergyDistribution; //TODO

            TH1F* ECalEnergyDistribution; //TODO
            TH1F* ECalBarrelEnergyDistribution;            
            TH1F* ECalEndcapEnergyDistribution;

            TH1F* ECalRecoEnergyDistribution; //TODO
            TH1F* ECalBarrelRecoEnergyDistribution; //TODO
            TH1F* ECalEndcapRecoEnergyDistribution; //TODO            
        
        //MCParticle
        TH1F* mcparticle_primaryparticleenergy;
        TH1F* mcparticle_nparticles;
        TH1F* mcparticle_nprimarys;
        TH1F* mcparticle_secondary_energies;
        TH1F* _MCParticleTotalEnergy;

        //HCalBarrel
        TH1F* hcalbarrel_nentries;
        TH1F* HCalBarrelFirstLayerOneHit;
        TH1F* HCalBarrelMomentaDistribution;        
        TProfile *HCalBarrelTotalEnergyPerLayer;
        TProfile *hcalbarrel_hits_vs_layer_profile;
        TProfile *hcalbarrel_allhits_vs_layer_profile;

        //HCalBarrelReco
        TH1F* HCalBarrelRecoNEntries; //TODO
        TH1F* HCalBarrelRecoFirstLayerTotalEnergy;
        TH1F* HCalBarrelRecoFirstLayerOneHit;        
        TProfile *HCalBarrelRecoTotalEnergyPerLayer;
        TProfile *hcalbarrelreco_hits_vs_layer_profile;
        //TProfile *hcalbarrelreco_allhits_vs_layer_profile;

        //HCalEndcap
        TH1F* hcalendcap_nentries;

        //HCalEndcapReco
        TH1F* HCalEndcapRecoNEntries; //TODO
        
        //ECalEndcap
        TH1F* ECalEndcapNEntries;

        //ECalEndcapReco
        TH1F* ECalEndcapRecoNEntries; //TODO
        
        //ECalBarrel
        TH1F* ECalBarrelNEntries;

        //ECalBarrelReco
        TH1F* ECalBarrelRecoNEntries; //TODO

       
        vector <string> colnames_to_analyze;

        int _nRun ;
        int _nEvt ;
};

#endif

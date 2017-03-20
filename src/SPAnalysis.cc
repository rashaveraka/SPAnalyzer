/*
 * author Ross McCoy, based on Christopher Milke's original instructional code
 * April 5, 2016
 */

/*
 *      TODO:
 *          Add another Marlin processor parameter that functions as multiplicative factor for energy distributions
 * 
 */

#include "SPAnalysis.h"
#include "IDDecoder.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/MCParticle.h>
#include <UTIL/CellIDDecoder.h>

#include <TFile.h>
#include <TH2D.h>

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"

#include <stdint.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
/*#Extra documentation for Ross.
#These processors all have more or less
#the same basic structure:
# 1) Constructor
#    This is where you register input parameters.
#    That is, you read in data from the steering xml.
#
# 2) init
#    Here you initialize any global variables,
#    and do any calculations needed BEFORE you
#    read in any of your slcio/stdhep data
#    files.
#
# 3) processEvent
#    This is the main part of the analysis code.
#    It takes as its argument an event from the 
#    input file you gave it. It will read one event
#    at a time, starting the function over again for
#    every event. Note that this means any values you
#    that you want to retain from one event to another
#    should be global to the file, as local variables
#    are reset for each event.
#
# 4) end
#    This is like init, but run after every event has
#    been processed. Use it to write out plot data or
#    calculate final statistics as needed.
#
/* There are a few other functions, like checkEvent,
   but I'm not terribly familiar with them, and have
   personally never needed to use them.*/


using namespace lcio;
using namespace marlin;
using namespace std;


SPAnalysis SPAnalysis;





SPAnalysis::SPAnalysis() : Processor("SPAnalysis") {
    // modify processor description
    _description = "Single Particle Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    
    MCParticleCollectionName="MCParticle";
    HCALBarrelCollectionName="HCalBarrelHits";
    HCALBarrelRecoCollectionName="HCalBarrelReco";
    HCALEndcapCollectionName="HCalEndcapHits";
    HCALEndcapRecoCollectionName="HCalEndcapReco";
    ECALBarrelCollectionName="ECalBarrelHits";
    ECALBarrelRecoCollectionName="ECalBarrelReco";
    ECALEndcapCollectionName="ECalEndcapHits";
    ECALEndcapRecoCollectionName="ECalEndcapReco";
    Detectorname = ""; 
    registerProcessorParameter("ROOTFileName","Name of the ROOT outputfile", ROOToutputfilename,string("output.root"));
    registerProcessorParameter("ParticleParameter","Particle name and energy", ParticleParameters,string("ERRNoInfo"));
    registerProcessorParameter("Multiplier","Multiplicative factor for histograms", Multiplier,float(1.0) );
}



void SPAnalysis::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    
    int ebins=1400;
    double emin=0;
    double emax=140;

    rootoutputfile=new TFile(ROOToutputfilename.c_str(),"RECREATE");
    if (rootoutputfile==NULL)
    {
    	cout<< "Can't open ROOT output: "<<ROOToutputfilename<<endl;
    }
    rootoutputfile->mkdir("TotalEnergyPlots");
    rootoutputfile->cd("TotalEnergyPlots");
    
        _MCParticleTotalEnergy = new TH1F("MCParticleTotalEnergy",std::string("Total Energy of all MCParticles for " + ParticleParameters).c_str(),300.0,0,30);

        HCalTotalEnergy = new TH1F("HCalTotalEnergy",std::string("Total Energy of all hits in the AHCal per event for " + ParticleParameters + ";Total HCal Energy per Event [GeV];# of Events at Binned Energy").c_str(),200,0,0.800*Multiplier);
        HCalRecoTotalEnergy = new TH1F("HCalRecoTotalEnergy",std::string("Total Energy of all hits in the AHCal per reco event for " + ParticleParameters + ";Total HCal Energy per Event [GeV];# of Events at Binned Energy").c_str(),200,0,25*Multiplier);
        HCalBarrelTotalEnergy = new TH1F("HCalBarrelTotalEnergy",std::string("Total Energy of all HCalBarrel Hits per event for " + ParticleParameters + ";Total HCalBarrel Energy per Event [GeV];# of Events at Binned Energy").c_str(),200,0,0.400*Multiplier);
        HCalEndcapTotalEnergy = new TH1F("HCalEndcapTotalEnergy",std::string("Total Energy of all HCalEndcap Hits per event for " + ParticleParameters + ";Total HCalEndcap Energy per Event [GeV];# of Events at Binned Energy").c_str(),200,0,0.100*Multiplier);
        HCalBarrelRecoTotalEnergy = new TH1F("HCalBarrelRecoTotalEnergy",std::string("Total Energy of all HCalBarrelReco Hits per event for " + ParticleParameters).c_str(),250,0,25*Multiplier);    
        HCalEndcapRecoTotalEnergy = new TH1F("HCalEndcapRecoTotalEnergy",std::string("Total Energy of all HCalEndcapReco Hits per event for " + ParticleParameters).c_str(),200,0,4*Multiplier);    

        ECalTotalEnergy = new TH1F("ECalTotalEnergy",std::string("Total Energy of all hits in the ECal per event for " + ParticleParameters + ";Total ECal Energy per Event [GeV];# of Events at Binned Energy").c_str(),200,0,0.800*Multiplier);
        ECalRecoTotalEnergy = new TH1F("ECalRecoTotalEnergy",std::string("Total Energy of all hits in the ECal per reco event for " + ParticleParameters + ";Total ECal Energy per Event [GeV];# of Events at Binned Energy").c_str(),200,0,15*Multiplier);
        ECalBarrelTotalEnergy = new TH1F("ECalBarrelTotalEnergy",std::string("Total energy of hits in ECalBarrel for " + ParticleParameters).c_str(),200,0,0.400);
        ECalEndcapTotalEnergy = new TH1F("ECalEndcapTotalEnergy",std::string("Total energy of hits in ECalEndcap for " + ParticleParameters).c_str(),200,0,0.400);    
        ECalBarrelRecoTotalEnergy = new TH1F("ECalBarrelRecoTotalEnergy",std::string("Total energy of hits in ECalBarrelReco per event for " + ParticleParameters).c_str(),200,0,2*Multiplier);
        ECalEndcapRecoTotalEnergy = new TH1F("ECalEndcapRecoTotalEnergy",std::string("Total energy of hits in ECalEndcapReco per event for " + ParticleParameters).c_str(),200,0,2*Multiplier);

        TotalEnergy = new TH1F("TotalEnergy",std::string("Total Energy of all hits in ECal and AHCal per event for " + ParticleParameters + ";Total Energy per Event [GeV];# of Events at Binned Energy").c_str(),200,0,0.6*Multiplier);//0.8 for 10 GeV Pions
        RecoTotalEnergy = new TH1F("RecoTotalEnergy",std::string("Total Energy of all hits in ECal and AHCal per reco event for " + ParticleParameters + ";Total Energy per Event [GeV];# of Events at Binned Energy").c_str(),200,0,20*Multiplier);//0.8 for 10 GeV Pions
        
        ECalvsHCalEnergy = new TH2F("HCalvsECalEnergy",std::string("Total Energy of all hits in ECal vs AHCal per reconstructed event for " + ParticleParameters + ";Total HCal Energy per Event [GeV];Total ECal Energy per Event [GeV]").c_str(),200,0,25*Multiplier,200,0,15*Multiplier);
    
    rootoutputfile->mkdir("EnergyDistributionPlots");
    rootoutputfile->cd("EnergyDistributionPlots");
        RecoTotalEnergyDistribution = new TH1F("RecoTotalEnergyDistribution",std::string("Energy Distribution of all Hits in ECal/AHCal per Post-Reco Event for " + ParticleParameters + ";Energy per Hit [GeV];# of Hits at Binned Energy").c_str(),200,0,0.3*Multiplier);//0.8 for 10 GeV Pions
        
        HCalBarrelEnergyDistribution = new TH1F("HCalBarrelEnergyDistribution",std::string("Distribution of deposited energy for each cell in HCalBarrel for " + ParticleParameters + ";Energy per Hit [GeV];# of Hits at Binned Energy").c_str(),1000,0,0.005);
        HCalEndcapEnergyDistribution = new TH1F("HCalEndcapEnergyDistribution",std::string("Distribution of deposited energy for each cell in HCalEndcap for " + ParticleParameters + ";Energy per Hit [GeV];# of Hits at Binned Energy").c_str(),1000,0,0.001);
        
        HCalBarrelRecoEnergyDistribution = new TH1F("HCalBarrelRecoEnergyDistribution",std::string("Distribution of deposited energy for each cell in HCalBarrelReco for " + ParticleParameters + ";Energy per Hit [GeV];# of Hits at Binned Energy").c_str(),500,0,1*Multiplier);        
        
        ECalEndcapEnergyDistribution = new TH1F("ECalEndcapEnergyDistribution",std::string("Distribution of deposited energy for each cell in ECalEndcap for " + ParticleParameters + ";Energy per Hit [GeV];# of Hits at Binned Energy").c_str(),1000,0,0.0005);
        ECalBarrelEnergyDistribution = new TH1F("ECalBarrelEnergyDistribution",std::string("Distribution of deposited energy for each cell in ECalBarrel for " + ParticleParameters + ";Energy per Hit [GeV];# of Hits at Binned Energy").c_str(),1000,0,0.001);
        
    rootoutputfile->mkdir("Calibration");
    rootoutputfile->cd("Calibration");
        DetectedMuonSignals = new TH1F("MuonCalibration", std::string("Distribution of detected muon signals in each cell in HCalBarrel for " + ParticleParameters).c_str(), 500, 0, 0.00500);
        DetectedMuonSignalsReco = new TH1F("MuonCalibrationReco", std::string("Distribution of detected muon signals in each cell in HCalBarrel for " + ParticleParameters + ";x-axis;y-axis").c_str(), 500, 0, 0.500);    

    rootoutputfile->mkdir("MCParticle");
    rootoutputfile->cd("MCParticle");
        mcparticle_nparticles  = new TH1F("MCParticleNParticle",std::string("MCParticleNParticle for " + ParticleParameters).c_str(),800, 0, 800 );
        mcparticle_primaryparticleenergy =new TH1F("MCParticlePrimaryParticleEnergy",std::string("MCParticlePrimaryParticleEnergy for " + ParticleParameters).c_str(), 2000, 9, 11 ); // set for muons right now
        mcparticle_nprimarys            	=new TH1F("MCParticleNPrimaryParticle",std::string("MCParticleNPrimaryParticle for " + ParticleParameters).c_str(), 20, 0, 20 );
        mcparticle_secondary_energies    =new TH1F("MCParticle Secondary Energies",std::string("MCParticle Secondary Energies for " + ParticleParameters).c_str(),800, 0, 0.04 );
    
    rootoutputfile->mkdir("HCALBarrel");
    rootoutputfile->cd("HCALBarrel");
        hcalbarrel_nentries	=new TH1F("HCALBarrelNentries",std::string("HCALBarrelNentries for " + ParticleParameters + ";# of HCalBarrel Hits per Event;# of Events with Binned Hits").c_str(),200,0,1000*Multiplier);
        hcalbarrel_hits_vs_layer_profile			= new TProfile("HCAL Barrel Hits_vs_Layer",std::string("HCAL Barrel Hits_vs_Layer for " + ParticleParameters).c_str(),HCAL_MAX_LAYERS,0,HCAL_MAX_LAYERS) ;
        hcalbarrel_allhits_vs_layer_profile		= new TProfile("HCAL Barrel AllHits_vs_Layer",std::string("HCAL Barrel AllHits_vs_Layer for " + ParticleParameters).c_str(),HCAL_MAX_LAYERS,0,HCAL_MAX_LAYERS) ;
        //HCalBarrelFirstLayerOneHit = new TH1F("HCalBarrelFirstLayerOneHit",std::string("Total Energy of Hits in HCalBarrel First Layer with Single Hit Detected for " + ParticleParameters, 600,0,6); // For 10GeV Pions
        HCalBarrelFirstLayerOneHit = new TH1F("HCalBarrelFirstLayerOneHit",std::string("Total Energy of Hits in First Layer of HCalBarrel with Single Hit Detected for " + ParticleParameters).c_str(), 500,0,0.05); // For 10GeV muon calibration
        HCalBarrelMomentaDistribution = new TH1F("HCalBarrelMomentaDistribution", std::string("Momenta distribution for each hit in HCalBarrel for " + ParticleParameters).c_str(),1000,0,1);
        HCalBarrelTotalEnergyPerLayer = new TProfile("HCalBarrelTotalEnergyPerLayer",std::string("Total energy deposited in HCalBarrel vs layers for " + ParticleParameters).c_str(),HCAL_MAX_LAYERS,0,HCAL_MAX_LAYERS);
    
    
    rootoutputfile->mkdir("HCALBarrelReco");
    rootoutputfile->cd("HCALBarrelReco");
        
        HCalBarrelRecoFirstLayerOneHit = new TH1F("HCalBarrelFirstLayerOneHit",std::string("Total Energy of Hits in First Layer of HCalBarrel with Single Hit Detected for " + ParticleParameters).c_str(), 600,0,3*Multiplier);    
        HCalBarrelRecoTotalEnergyPerLayer = new TProfile("HCalBarrelRecoTotalEnergyPerLayer",std::string("Total energy deposited in HCalBarrelReco vs layers for " + ParticleParameters).c_str(),HCAL_MAX_LAYERS,0,HCAL_MAX_LAYERS);
        hcalbarrelreco_hits_vs_layer_profile = new TProfile("HCALBarrelRecoHits_vs_Layer",std::string("HCAL Barrel Reco Hits_vs_Layer for " + ParticleParameters).c_str(),HCAL_MAX_LAYERS,0,HCAL_MAX_LAYERS);
        //hcalbarrelreco_allhits_vs_layer_profile = new TProfile("HCALBarrelRecoAllHits_vs_Layer",std::string("HCAL Barrel Reco AllHits_vs_Layer for " + ParticleParameters).c_str(),HCAL_MAX_LAYERS,0,HCAL_MAX_LAYERS);
        HCalBarrelRecoFirstLayerTotalEnergy = new TH1F("HCalBarrelRecoFirstLayerTotalEnergy",std::string("Total Energy of Hits in First Layer of HCalBarrelReco per Event for " + ParticleParameters).c_str(),200,0,6*Multiplier);
        //HCalBarrelRecoMomentaDistribution = new TH1F("HCalBarrelRecoMomentaDistribution", std::string("Momenta distribution for each hit in HCalBarrelReco for " + ParticleParameters).c_str(),1000,0,1);
    
    
    rootoutputfile->mkdir("HCALEndcap");
    rootoutputfile->cd("HCALEndcap");
        hcalendcap_nentries	=new TH1F("HCALEndcapNentries",std::string("HCALEndcapNentries for " + ParticleParameters).c_str(),200,0,5000);

    rootoutputfile->mkdir("HCALEndcapReco");
    rootoutputfile->cd("HCALEndcapReco");    
         //TODO: Add relevant histograms   
    
    rootoutputfile->mkdir("ECALEndcap");
    rootoutputfile->cd("ECALEndcap");
        ECalEndcapNEntries	=new TH1F("ECALEndcapNentries",std::string("ECALEndcapNentries for " + ParticleParameters).c_str(),200,0,5000);

    rootoutputfile->mkdir("ECALEndcapReco");
    rootoutputfile->cd("ECALEndcapReco");    
        //TODO: Add relevant histograms
    
    rootoutputfile->mkdir("ECALBarrel");
    rootoutputfile->cd("ECALBarrel");
        ECalBarrelNEntries	=new TH1F("ECALBarrelNentries",std::string("ECALBarrelNentries for " + ParticleParameters).c_str(),200,0,5000);

    rootoutputfile->mkdir("ECALBarrelReco");
    rootoutputfile->cd("ECALBarrelReco");
        //TODO: Add relevant histograms
    
   
    colnames_to_analyze.push_back(MCParticleCollectionName);
    colnames_to_analyze.push_back(HCALBarrelCollectionName);
    colnames_to_analyze.push_back(HCALBarrelRecoCollectionName);
    colnames_to_analyze.push_back(HCALEndcapCollectionName);
    colnames_to_analyze.push_back(HCALEndcapRecoCollectionName);    
    colnames_to_analyze.push_back(ECALBarrelCollectionName);
    colnames_to_analyze.push_back(ECALBarrelRecoCollectionName);    
    colnames_to_analyze.push_back(ECALEndcapCollectionName);
    colnames_to_analyze.push_back(ECALEndcapRecoCollectionName);    

    printParameters();

    _nRun = 0 ;
    _nEvt = 0 ;
  
}



void SPAnalysis::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void SPAnalysis::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...
    /* Commenting out for development
    LCCollection* col = evt->getCollection( _colName ) ;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
        int cellID = 0;
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(hitIndex) );
           _hitmap->Fill(hit->getEnergy());
	   cellID = hit->getCellID0();
	   cout << " Cell ID is: " << cellID << endl;
        } 
        _hitmap->Fill(hit->getEnergy());
    }*/
    
    const vector <string>*  colNames =  evt->getCollectionNames(); // all collection names in the file
    LCCollection* col = NULL;
    bool has_hcal_barrel_hits=false;
    float HCalTotalEnergySum = 0.0;
    float HCalRecoTotalEnergySum = 0.0;
    float ECalTotalEnergySum = 0.0;
    float ECalRecoTotalEnergySum = 0.0;
    float TotalEnergySum = 0.0;
    float RecoTotalEnergySum = 0.0;
    
    int cellhitarray[HCAL_MAX_LAYERS];
    int totalhitarray[HCAL_MAX_LAYERS];
    float TotalEnergyLayers[HCAL_MAX_LAYERS];
    
    
    string encodingString;

    
    unsigned int n_mc_primaries=0;
    
    for( vector <string>::const_iterator it = colnames_to_analyze.begin(); it != colnames_to_analyze.end() ; it++ )
    {
        try{
                col = evt->getCollection( *it);
        }
        catch( lcio::DataNotAvailableException & e )
        {
                streamlog_out(WARNING) << *it << " collection not available" << std::endl;
                col = NULL;
        }
        
        // this will only be entered if the collection is available
        if( col != NULL )
        {
            int entries = col->getNumberOfElements();
            streamlog_out(DEBUG) <<*it<< " "<<col->getTypeName()<<" "<<entries<<endl;

                if ( ((*it).compare(MCParticleCollectionName)==0)&& (col->getTypeName()=="MCParticle") )
                {
                    mcparticle_nparticles->Fill(entries);

                    for (int i=0; i< entries; i++)
                    {
                        MCParticle* myparticle=  dynamic_cast<MCParticle*>( col->getElementAt(i));

                        if ((myparticle->getParents().size()==0) &&(myparticle->getDaughters().size()>0))//thats a primary
                        {
                            mcparticle_primaryparticleenergy->Fill(myparticle->getEnergy());


                            n_mc_primaries++;

                        }
                        else //look at the secondaries
                        {
                            mcparticle_secondary_energies->Fill(myparticle->getEnergy());
                        }

                    }
                    mcparticle_nprimarys->Fill(n_mc_primaries);
                } //end of MCParticle block                
        
    
   		if ( ((*it).compare(HCALBarrelCollectionName)==0)&& (col->getTypeName()=="SimCalorimeterHit") ) //HCAL Barrel
    		{
                        encodingString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
                        streamlog_out(MESSAGE) << "Encoding string for HCalBarrel " << encodingString << endl;
                        
                        int numHitsHCalBarrelFirstLayer = 0;
                        float HCalBarrelFirstLayerTotalEnergySum = 0.0;
                        float HCalBarrelTotalEnergySum = 0.0;
                        
                        UTIL::BitField64 cellIDDecoder(encodingString);

    			if (entries>0)
    			{

    				hcalbarrel_nentries->Fill(entries);
    				has_hcal_barrel_hits=true;

    				for (int i=0;i<HCAL_MAX_LAYERS; i++)
    				{
    				    	cellhitarray[i]=0;
    				    	totalhitarray[i]=0;
                                        TotalEnergyLayers[i]=0;
    				}
                                
    				for (int i=0; i< entries; i++)
    				{
    				    	SimCalorimeterHit * myhit=  dynamic_cast<SimCalorimeterHit*>( col->getElementAt(i));
                                        HCalBarrelTotalEnergySum += myhit->getEnergy(); // Need to check and see if above 0.6MeV?
                                        //HCalBarrelMomentaDistribution->Fill(myhit->getMomentum());
                                        
                                        DD4hep::long64 id = myhit->getCellID0();
                                        // where *it points to the hit
                                        cellIDDecoder.setValue( id );
                                        int layer = cellIDDecoder["layer"].value();
					//streamlog_out(MESSAGE) << " Hit in HCalBarrel Layer " << layer << " with ID " << id << endl;
                                        
                                        if(layer == 1)
                                        {
                                            numHitsHCalBarrelFirstLayer++;
                                            HCalBarrelFirstLayerTotalEnergySum += myhit->getEnergy();
                                        }
                                        
                                        cellhitarray[layer]++;
                                        TotalEnergyLayers[layer]+=myhit->getEnergy();
                                        
                                        
    				        // Disabled for now, used for muon calibration
                                        DetectedMuonSignals->Fill(myhit->getEnergy());
                                        HCalBarrelEnergyDistribution->Fill(myhit->getEnergy()); //Need to check and see if above 0.6MeV?
                                
                                }
                                HCalTotalEnergySum += HCalBarrelTotalEnergySum;
                                HCalBarrelTotalEnergy->Fill(HCalBarrelTotalEnergySum);
                                TotalEnergySum += HCalBarrelTotalEnergySum;
                                
                                if(numHitsHCalBarrelFirstLayer == 1)
                                {
                                    HCalBarrelFirstLayerOneHit->Fill(HCalBarrelTotalEnergySum);
                                }
                                
                                for (int i=0;i<HCAL_MAX_LAYERS; i++)
    				{
                                    HCalBarrelTotalEnergyPerLayer->Fill(i,TotalEnergyLayers[i]);
                                    hcalbarrel_hits_vs_layer_profile->Fill(i,cellhitarray[i]);
    				}

    			}

    		}
                if ( ((*it).compare(HCALBarrelRecoCollectionName)==0)&& (col->getTypeName()=="CalorimeterHit") ) //HCAL Barrel Reco
    		{
                        encodingString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
                        streamlog_out(MESSAGE) << "Encoding string for HCalBarrelReco " << encodingString << endl;
                        
                        int numHitsHCalBarrelRecoFirstLayer = 0;
                        float HCalBarrelRecoFirstLayerTotalEnergySum = 0.0;
                        float HCalBarrelRecoTotalEnergySum = 0.0;                        
                        
                        UTIL::BitField64 cellIDDecoder(encodingString);
                        
    			if (entries>0)
    			{
    				for (int i=0;i<HCAL_MAX_LAYERS; i++)
    				{
    				    	cellhitarray[i]=0;
    				    	totalhitarray[i]=0;
                                        TotalEnergyLayers[i]=0;
    				}

    				for (int i=0; i< entries; i++)
    				{
    				    	CalorimeterHit * myhit=  dynamic_cast<CalorimeterHit*>( col->getElementAt(i));
                                        DD4hep::long64 id = myhit->getCellID0();
                                        cellIDDecoder.setValue( id );
                                        int layer = cellIDDecoder["layer"].value();
                                        if(layer == 1)
                                        {
                                            numHitsHCalBarrelRecoFirstLayer++;
                                            HCalBarrelRecoFirstLayerTotalEnergySum += myhit->getEnergy();
                                        }
                                        //HCalBarrelRecoMomentaDistribution->Fill(myhit->getMomentum()); No getMomentum() function in this collection
                                        TotalEnergyLayers[layer] += myhit->getEnergy();
                                        cellhitarray[layer]++;
   				    	HCalBarrelRecoTotalEnergySum += myhit->getEnergy();

    				        // Disabled for now, used for muon calibration
                                        //DetectedMuonSignalsReco->Fill(myhit->getEnergy());
                                        HCalBarrelRecoEnergyDistribution->Fill(myhit->getEnergy()); //Need to check and see if above 0.6MeV?
                                        RecoTotalEnergyDistribution->Fill(myhit->getEnergy());
                                }
                                HCalBarrelRecoFirstLayerTotalEnergy->Fill(HCalBarrelRecoFirstLayerTotalEnergySum);
                                HCalBarrelRecoTotalEnergy->Fill(HCalBarrelRecoTotalEnergySum);
                                HCalRecoTotalEnergySum += HCalBarrelRecoTotalEnergySum;
                                RecoTotalEnergySum += HCalBarrelRecoTotalEnergySum;                               
                                
                                if(numHitsHCalBarrelRecoFirstLayer == 1)
                                {
                                    HCalBarrelRecoFirstLayerOneHit->Fill(HCalBarrelRecoFirstLayerTotalEnergySum);
                                }                                
                                
                                for (int i=0;i<HCAL_MAX_LAYERS; i++)
    				{
                                    HCalBarrelRecoTotalEnergyPerLayer->Fill(i,TotalEnergyLayers[i]);
                                    hcalbarrelreco_hits_vs_layer_profile->Fill(i,cellhitarray[i]);
    				}

    			}

    		}

            
            
   		if ( ((*it).compare(HCALEndcapCollectionName)==0)&& (col->getTypeName()=="SimCalorimeterHit") ) //HCAL Endcap
    		{
    			if (entries>0)
    			{

    				hcalendcap_nentries->Fill(entries);
                                float HCalEndcapTotalEnergySum = 0.0;

    				for (int i=0; i< entries; i++)
    				{
    				    	SimCalorimeterHit * myhit=  dynamic_cast<SimCalorimeterHit*>( col->getElementAt(i));
    				    	
                                        HCalEndcapTotalEnergySum += myhit->getEnergy(); // Need to check and see if above 0.6MeV?
                                        

    				       // Disabled for now, used for muon calibration
                                       // DetectedMuonSignals->Fill(myhit->getEnergy());
                                        HCalEndcapEnergyDistribution->Fill(myhit->getEnergy()); //Need to check and see if above 0.6MeV?
                                
                                }
                                HCalEndcapTotalEnergy->Fill(HCalEndcapTotalEnergySum);
                                HCalTotalEnergySum += HCalEndcapTotalEnergySum;
                                TotalEnergySum += HCalEndcapTotalEnergySum;

    			}

    		}

   		if ( ((*it).compare(HCALEndcapRecoCollectionName)==0)&& (col->getTypeName()=="CalorimeterHit") ) //HCAL Endcap
    		{
    			if (entries>0)
    			{

    				//hcalendcap_nentries->Fill(entries);
                                float HCalEndcapRecoTotalEnergySum = 0.0;

    				for (int i=0; i< entries; i++)
    				{
    				    	CalorimeterHit * myhit=  dynamic_cast<CalorimeterHit*>( col->getElementAt(i));
    				    	
                                        HCalEndcapRecoTotalEnergySum += myhit->getEnergy(); // Need to check and see if above 0.6MeV?
                                        

    				       // Disabled for now, used for muon calibration
                                       // DetectedMuonSignals->Fill(myhit->getEnergy());
                                       // HCalEndcapEnergyDistribution->Fill(myhit->getEnergy()); //Need to check and see if above 0.6MeV?
                                        RecoTotalEnergyDistribution->Fill(myhit->getEnergy());
                                }
                                HCalEndcapRecoTotalEnergy->Fill(HCalEndcapRecoTotalEnergySum);
                                HCalRecoTotalEnergySum += HCalEndcapRecoTotalEnergySum;
                                RecoTotalEnergySum += HCalEndcapRecoTotalEnergySum;

    			}

    		}
            
   		if ( ((*it).compare(ECALEndcapCollectionName)==0)&& (col->getTypeName()=="SimCalorimeterHit") ) //ECAL Endcap
    		{
    			if (entries>0)
    			{

    				ECalEndcapNEntries->Fill(entries);
                                float ECalEndcapTotalEnergySum = 0.0;

    				for (int i=0; i< entries; i++)
    				{
    				    	SimCalorimeterHit * myhit=  dynamic_cast<SimCalorimeterHit*>( col->getElementAt(i));
    				    	
                                        ECalEndcapTotalEnergySum += myhit->getEnergy(); // Need to check and see if above 0.6MeV?
                                        

    				       // Disabled for now, used for muon calibration
                                       // DetectedMuonSignals->Fill(myhit->getEnergy());
                                        ECalEndcapEnergyDistribution->Fill(myhit->getEnergy()); //Need to check and see if above 0.6MeV?
                                
                                }
                                ECalEndcapTotalEnergy->Fill(ECalEndcapTotalEnergySum);
                                ECalTotalEnergySum += ECalEndcapTotalEnergySum;
                                TotalEnergySum += ECalEndcapTotalEnergySum;

    			}

    		}
            
   		if ( ((*it).compare(ECALEndcapRecoCollectionName)==0)&& (col->getTypeName()=="CalorimeterHit") ) //ECAL Endcap
    		{
    			if (entries>0)
    			{

    				//ECalEndcapNEntries->Fill(entries);
                                float ECalEndcapRecoTotalEnergySum = 0.0;

    				for (int i=0; i< entries; i++)
    				{
    				    	CalorimeterHit * myhit=  dynamic_cast<CalorimeterHit*>( col->getElementAt(i));
    				    	
                                        ECalEndcapRecoTotalEnergySum += myhit->getEnergy(); // Need to check and see if above 0.6MeV?
                                        

    				       // Disabled for now, used for muon calibration
                                       // DetectedMuonSignals->Fill(myhit->getEnergy());
                                       //ECalEndcapEnergyDistribution->Fill(myhit->getEnergy()); //Need to check and see if above 0.6MeV?
                                        RecoTotalEnergyDistribution->Fill(myhit->getEnergy());
                                }
                                ECalEndcapRecoTotalEnergy->Fill(ECalEndcapRecoTotalEnergySum);
                                ECalRecoTotalEnergySum += ECalEndcapRecoTotalEnergySum;
                                RecoTotalEnergySum += ECalEndcapRecoTotalEnergySum;

    			}

    		}            
            
   		if ( ((*it).compare(ECALBarrelCollectionName)==0)&& (col->getTypeName()=="SimCalorimeterHit") ) //ECAL Barrel
    		{
    			if (entries>0)
    			{

    				ECalBarrelNEntries->Fill(entries);
                                float ECalBarrelTotalEnergySum = 0.0;

    				for (int i=0; i< entries; i++)
    				{
    				    	SimCalorimeterHit * myhit=  dynamic_cast<SimCalorimeterHit*>( col->getElementAt(i));
    				    	
                                        ECalBarrelTotalEnergySum += myhit->getEnergy(); // Need to check and see if above 0.6MeV?
                                        

    				       // Disabled for now, used for muon calibration
                                       // DetectedMuonSignals->Fill(myhit->getEnergy());
                                        ECalBarrelEnergyDistribution->Fill(myhit->getEnergy()); //Need to check and see if above 0.6MeV?
                                
                                }
                                ECalBarrelTotalEnergy->Fill(ECalBarrelTotalEnergySum);
                                ECalTotalEnergySum += ECalBarrelTotalEnergySum;
                                TotalEnergySum += ECalBarrelTotalEnergySum;

    			}

    		}

   		if ( ((*it).compare(ECALBarrelRecoCollectionName)==0)&& (col->getTypeName()=="CalorimeterHit") ) //ECAL Barrel
    		{
    			if (entries>0)
    			{

                                float ECalBarrelRecoTotalEnergySum = 0.0;

    				for (int i=0; i< entries; i++)
    				{
    				    	CalorimeterHit * myhit=  dynamic_cast<CalorimeterHit*>( col->getElementAt(i));
    				    	
                                        ECalBarrelRecoTotalEnergySum += myhit->getEnergy(); // Need to check and see if above 0.6MeV?
                                        

    				       // Disabled for now, used for muon calibration
                                       // DetectedMuonSignals->Fill(myhit->getEnergy());
                                       // ECalBarrelEnergyDistribution->Fill(myhit->getEnergy()); //Need to check and see if above 0.6MeV?
                                        RecoTotalEnergyDistribution->Fill(myhit->getEnergy());
                                }
                                ECalBarrelRecoTotalEnergy->Fill(ECalBarrelRecoTotalEnergySum);
                                ECalRecoTotalEnergySum += ECalBarrelRecoTotalEnergySum;
                                RecoTotalEnergySum += ECalBarrelRecoTotalEnergySum;

    			}

    		}            
    
        }
    }
    
        ECalvsHCalEnergy->Fill(HCalRecoTotalEnergySum,ECalRecoTotalEnergySum);
        HCalTotalEnergy->Fill(HCalTotalEnergySum);
        HCalRecoTotalEnergy->Fill(HCalRecoTotalEnergySum);
        ECalTotalEnergy->Fill(ECalTotalEnergySum);
        ECalRecoTotalEnergy->Fill(ECalRecoTotalEnergySum);
        TotalEnergy->Fill(TotalEnergySum);
        RecoTotalEnergy->Fill(RecoTotalEnergySum);
        streamlog_out(MESSAGE) << "   processing event: " << evt->getEventNumber()
        << "   in run:  " << evt->getRunNumber() << std::endl ;
        
    _nEvt ++ ;
}



void SPAnalysis::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void SPAnalysis::end(){ 
    rootoutputfile->Write();
}

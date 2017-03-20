print "Starting Analysis.  This may take a moment..."


from pyLCIO import IOIMPL, EVENT,UTIL
from pyLCIO.EVENT import DataNotAvailableException
from pyLCIO.io.LcioReader import LcioReader
from array import array
from ROOT import *


# Input File *change to slcio file that will analyized 
slcioFile = '5GeV_65_pions_10000events_reco.slcio'

# output root file *change sample.root to desired name
rootFile = TFile("sample.root","recreate")

# create Histograms
HitEnergyHistogram = TH1D( 'TotalDepositedEnergy', 'Event Energy Deposit;Hit Energy [GeV];Entries', 100, 0., 1.)
HitsvsLayerHistogram = TH1D("HitsvsLayers","Hits vs. Layers;Layers;HCAL Barrel Hits", 40, 0., 39.)
NumberOfHitsPerEventHistogram = TH1D( 'HitsPerEvent', 'Hits Per Event;Hits;Events', 400, 0., 1600.)
HitEnergyHistogram_Digi_Reco = TH1D( 'TotalDepositedEnergyReco', 'Event Energy Deposit DigiReco;Hit Energy [GeV];Entries', 100, 0., 10.)


# create a reader and open the slcio file)
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(slcioFile)

# loop over all events in the file
for event in reader:
    
    print "Collecting Hits in HCAL Barrel, HCAL End Cap, ECAL Barrel, ECAL End Cap and summing energies for event %s" % ( event.getEventNumber() )
    
    hitEnergyTotal = 0
    numberofHits = 0
    hitEnergyTotal_Digi_Reco = 0

    #Get individual hit collections
    HCALBarrelhitCollection = event.getCollection( 'HCalBarrelHits' )
    HCALEndCaphitCollection = event.getCollection( 'HCalEndcapHits' )
    ECALBarrelhitCollection = event.getCollection( 'ECalBarrelHits' )
    ECALEndCaphitCollection = event.getCollection( 'ECalEndcapHits' )
    
    # get the cell ID encoding string from the collection parameters
    cellIdEncoding_HCALBarrel = HCALBarrelhitCollection.getParameters().getStringVal( EVENT.LCIO.CellIDEncoding )
    cellIdEncoding_HCALEndCap = HCALEndCaphitCollection.getParameters().getStringVal( EVENT.LCIO.CellIDEncoding )
    cellIdEncoding_ECALBarrel = ECALBarrelhitCollection.getParameters().getStringVal( EVENT.LCIO.CellIDEncoding )
    cellIdEncoding_ECALEndCap = ECALEndCaphitCollection.getParameters().getStringVal( EVENT.LCIO.CellIDEncoding )
    
    # define a cell ID decoder for the collection
    idDecoder_HCALBarrel = UTIL.BitField64( cellIdEncoding_HCALBarrel )
    idDecoder_HCALEndCap = UTIL.BitField64( cellIdEncoding_HCALEndCap )
    idDecoder_ECALBarrel = UTIL.BitField64( cellIdEncoding_ECALBarrel )
    idDecoder_ECALEndCap = UTIL.BitField64( cellIdEncoding_ECALEndCap )
    
    
    # Cycle through each hit in the collections and record the amount of hits and sum their energies of the event 
    for hit in HCALBarrelhitCollection:
        # combine the two 32 bit cell IDs of the hit into one 64 bit integer
        cellID = long( hit.getCellID0() & 0xffffffff ) | ( long( hit.getCellID1() ) << 32 )
        # set up the ID decoder for this cell ID
        idDecoder_HCALBarrel.setValue( cellID )
        # access the field information using a valid field from the cell ID encoding string and fill histogram with layer info
        HitsvsLayerHistogram.Fill(idDecoder_HCALBarrel['layer'].value())
        
        numberofHits += 1
        hitEnergyTotal+=hit.getEnergy()

    for hit in HCALEndCaphitCollection:
        # combine the two 32 bit cell IDs of the hit into one 64 bit integer
        cellID = long( hit.getCellID0() & 0xffffffff ) | ( long( hit.getCellID1() ) << 32 )
        # set up the ID decoder for this cell ID
        idDecoder_HCALEndCap.setValue( cellID )
        # access the field information using a valid field from the cell ID encoding string and fill histogram with layer info
        HitsvsLayerHistogram.Fill(idDecoder_HCALEndCap['layer'].value())
        
        numberofHits += 1
        hitEnergyTotal+=hit.getEnergy()

    for hit in ECALBarrelhitCollection:
        # combine the two 32 bit cell IDs of the hit into one 64 bit integer
        cellID = long( hit.getCellID0() & 0xffffffff ) | ( long( hit.getCellID1() ) << 32 )
        # set up the ID decoder for this cell ID
        idDecoder_ECALBarrel.setValue( cellID )
        # access the field information using a valid field from the cell ID encoding string and fill histogram with layer info
        HitsvsLayerHistogram.Fill(idDecoder_ECALBarrel['layer'].value())
        
        numberofHits += 1
        hitEnergyTotal+=hit.getEnergy()

    for hit in ECALEndCaphitCollection:
        # combine the two 32 bit cell IDs of the hit into one 64 bit integer
        cellID = long( hit.getCellID0() & 0xffffffff ) | ( long( hit.getCellID1() ) << 32 )
        # set up the ID decoder for this cell ID
        idDecoder_ECALEndCap.setValue( cellID )
        # access the field information using a valid field from the cell ID encoding string and fill histogram with layer info
        HitsvsLayerHistogram.Fill(idDecoder_ECALEndCap['layer'].value())
        
        numberofHits += 1
        hitEnergyTotal+=hit.getEnergy()

    # Cycle through each hit in the digi/reco collections and sum their energies
	if "ECalBarrelReco" in event.getCollectionNames(): 
		hits = event.getCollection("ECalEndcapReco")
		for hit in hits:
			hitEnergyTotal_Digi_Reco+=hit.getEnergy()
	if "ECalEndcapReco" in event.getCollectionNames(): 
		hits = event.getCollection("ECalEndcapReco")
		for hit in hits:
			hitEnergyTotal_Digi_Reco+=hit.getEnergy()
	if "HCalEndcapReco" in event.getCollectionNames(): 
		hits = event.getCollection("HCalEndcapReco")
		for hit in hits:
			hitEnergyTotal_Digi_Reco+=hit.getEnergy()		
	if "HCalBarrelReco" in event.getCollectionNames(): 
		hits = event.getCollection("HCalBarrelReco")
		for hit in hits:
			hitEnergyTotal_Digi_Reco+=hit.getEnergy()
        

    # Fill the histograms 
    HitEnergyHistogram.Fill(hitEnergyTotal)
    NumberOfHitsPerEventHistogram.Fill(numberofHits)
    HitEnergyHistogram_Digi_Reco.Fill(hitEnergyTotal_Digi_Reco)
reader.close()
    
# Fit histograms with Gaussian       ~ working on errors~
#HitEnergyHistogram.Fit("gaus")

# write and close the root file
rootFile.Write()
rootFile.Close()

print "Root file created: "+str(rootFile)





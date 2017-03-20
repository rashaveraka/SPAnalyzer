/*
 * IDDecoder.cpp
 *
 *  Created on: Jul 10, 2014
 *      Author: stanitz
 */

#include "IDDecoder.h"
#include "marlin/Processor.h"

#include <iostream>


IDDecoder::IDDecoder(std::string detectorname,std::string detectortype)
{
	Layer=0;
	internal_id=0;
	subdetectorname_id=0;
	detectorname_id=1; //only one detector for our purposes right now




	if (detectortype.compare("HCALBarrel")==0)
	{
			subdetectorname_id=4;
	}
    if (detectortype.compare("ECALBarrel")==0)
	{
			subdetectorname_id=3;
	}
    if (detectortype.compare("MuonBarrel")==0)
 	{
 			subdetectorname_id=5;
 	}


	internal_id=1000*detectorname_id+subdetectorname_id;
	streamlog_out(MESSAGE) <<detectorname<< " Internal ID " <<detectorname_id<<" "<<internal_id<<std::endl;

}

IDDecoder::~IDDecoder()
{
	// TODO Auto-generated destructor stub
}

void IDDecoder::SetCellID(int id0, int id1)
{

	//HCAL Barrel
    streamlog_out(MESSAGE) << "SetCellID called with ID0 " << id0 << " ID1 " << id1 << std::endl;
	if (internal_id==1004) // hcal barrel
	{
		Layer= (id0 >> 13)&0x3f;
                streamlog_out(MESSAGE) << "Layer is " << Layer << std::endl;
	}
//Muon Barrel
	if (internal_id==1005) // muon barrel
	{
		Layer= (id0 >> 15)&0xff;
	}

//ECAL Barrel

	if (internal_id==1003) // ecal barrel
	{
		Layer= (id0 >> 13)&0x3f;
	}



}



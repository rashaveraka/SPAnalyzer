/*
 * IDDecoder.h
 *
 *  Created on: Jul 10, 2014
 *      Author: stanitz
 */

#ifndef IDDECODER_H_
#define IDDECODER_H_

#include <string>

class IDDecoder
{
public:
	IDDecoder(std::string detectorname,std::string detectortype);
	virtual ~IDDecoder();
	void SetCellID(int id0, int id1);

	int getBarrel() const
	{
		return barrel;
	}

	int getDetectorId() const
	{
		return DetectorID;
	}

	int getLayer() const
	{
		return Layer;
	}

private:
	int Layer;
	int DetectorID;
	int barrel;
	int id0,id1;
	std::string detectorname;
	int detectorname_id;
	int subdetectorname_id; //1: Vertex: 2 : Tracker 3: ECAL 4: HCAL 5: Muons / x10 10 Endcap
	int internal_id;
};

#endif /* IDDECODER_H_ */

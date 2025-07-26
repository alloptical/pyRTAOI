// LibraryName parameters
#ifndef PARAMETER_H
#define PARAMETER_H

#include <GenerateHologramCUDA.h>
#include <vector>
#include "NIDAQmx.h"

namespace pa{
	class parameters;
}
class parameters{
public:
	std::string file_name;
	unsigned char *h_pSLM;
	int N_iterations, method, num_spots, PowerPerCell, photoDuration, photoSamples;
	float* x, *y, *z, *I, *h_obtainedAmps, *h_test,*x_i,*y_i;
	unsigned char *phaseMaskToSLM;
	char* insert_text;
	bool FlagInvert;
	double Ta[10], Tb[10];
	float64 photoVolt, daqSampleRate;


	// USING VECTORS:
	std::vector<float> xx, yy, zz, II, xx_i, yy_i;

	parameters(){
		h_pSLM = new unsigned char[512 * 512];
		N_iterations = 30;
		h_obtainedAmps = new float[512 * 512];
		//h_test = new float[512 * 512]; // added 20180424
		phaseMaskToSLM = new unsigned char[512 * 512];
		FlagInvert = false;
		method = 2;
		PowerPerCell = 6;
		photoVolt = 0;
		photoDuration = 10;// in ms
		daqSampleRate = 1000.0;
		photoSamples = 101; 

	};

	~parameters(){
		delete h_obtainedAmps;
		delete phaseMaskToSLM;
		delete h_pSLM;
		//delete h_test;
	};

	friend class holoblink;
};

#endif PARAMETER_H
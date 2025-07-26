
#include "holoblink.h"
#include <QApplication>





int main(int argc, char *argv[])
{
	// for heap corruption diagnose
	_CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF);


	// parameers for computing hologram are defined in parameter.h constructor

    //start cuda
	float *h_pSLMstart = new float[512*512];
	h_pSLMstart = { 0 };
	int cuda_error1 = startCUDA(h_pSLMstart, 0);
	delete h_pSLMstart;

	QApplication a(argc, argv);
    HoloBlink w;
    w.show();

    return a.exec();
}

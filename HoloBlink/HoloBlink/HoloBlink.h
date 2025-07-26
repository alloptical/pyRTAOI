#ifndef HOLOBLINK_H
#define HOLOBLINK_H

#include <QMainWindow>
#include <GenerateHologramCUDA.h>
#include <parameters.h>
#include <vector>
#include <QTableWidget>
#include <opencv2/core/core.hpp>>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <QTcpSocket>
#include <QTcpServer>
#include <QDebug>
#include "Blink_SDK.H" 
#include "Logger.h"

namespace Ui {
class HoloBlink;
}

class HoloBlink : public QMainWindow
{
    Q_OBJECT

public:

	parameters par;
	QTcpServer *server;
	QTcpSocket *client;
	//daq writing task  
	TaskHandle  WriteTaskHandle;

	int32 written;
	/*DAQ config*/
	std::vector<float64> WriteData;

	// logger
	CPlusPlusLogging::Logger* pLogger; // Create the object pointer for Logger Class

	int HoloBlink::invertInputs(float*, int);


    explicit HoloBlink(QWidget *parent = 0);
    ~HoloBlink();

	int get_hologram();
	void updateIinTable();
	float getTotalPowerOnSample();
	void transformXY(float*, float*, double*, double*, int);
	void transformXY_vect(); // added 20180424
	int readMessage(char, int);
	void updateSpots();
	void updateImages();
	cv::Mat make_save_img(float* , float*, bool );
	void weightIntensity();

private slots:

    void on_actionOpen_triggered();
    void on_actionSave_triggered();
    void on_transform_pushButton_clicked();
	int on_insert_pushButton_clicked();
	void on_invert_radioButton_clicked();
	void on_powerPerCell_valueChanged();
	void on_putOnSLM_pushButton_clicked();
	void on_InitialiseSLM_pushButton_clicked();

	void acceptConnection();
	void startRead();


private:
    Ui::HoloBlink *ui;
	QTableWidget* m_pTableWidget;
	QStringList m_TableHeader;
	static const QColor RED;
	QPixmap  BLANKIMG;

};

#endif // HOLOBLINK_H

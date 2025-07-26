/*
TO DO:
 - clean up for loading centroid bit (not using float* for xyzI anymore)
    0. deal with weighting
	1.show targets before transform
	2.save all target and phase masks during experiment
		
*/




#include "holoblink.h"
#include "ui_holoblink.h"
#include <QFileDialog>
#include <QFile>
#include <QMessageBox>
#include <QTextStream>
#include <QGraphicsItem>
#include <iostream>
#include <string.h>
#include <QtGui>
#include "mat.h"
#include "matrix.h"
#include "mex.h"
#include <exception>
#include <ctime>



using namespace cv;
using namespace CPlusPlusLogging; // Logger for error log

// global vars
std::vector<float> x_spots;
std::vector<float> y_spots;
std::vector<float> I;
// weighting -- copied from weighting file in settings.yml
double p_weight[2] = { -0.635949168274612, 316.724831762549};
double SteepnessFudgeFactor = 2;

// transform matrix -- import from file 
double Ta[3];
double Tb[3];
double Tc[3];
float XOFFSET = 5; // xoffset>0 moves spots to the right - this is opposite to the value in yml file
float YOFFSET = 0; // same as in yml YOFFSET>0 moves them down
bool FLAG_XY_CHANGE = FALSE;
bool FLAG_ICHANGE = FALSE;
bool FLAG_TRIGGER_PHOTO = FALSE;
bool FLAG_CONNECTED = FALSE;
bool FLAG_DURATION_UPDATED = FALSE;

// change this file name for new transform
const char *transform_file = "C:/Users/zoo/Dropbox/Bruker1/SLM/Transforms/ZZ_2DTransform_2020924_nikon16x1p14x_OLD_T.mat";
const char *transform_name = "tform2T";

// for test only
//const char *transform_file = "C:/Data/Zoe/VSprojects/HoloBlink/test files/test_tform.mat";
//const char *transform_name = "ans";

// Blink Configure Code
/*********************************************/
const int          board_number = 1;
const unsigned int bits_per_pixel = 8U;
const unsigned int pixel_dimension = 512U;
const bool         is_nematic_type = true;
const bool         RAM_write_enable = true;
const bool         use_GPU_if_available = true;
const char* const  regional_lut_file = "C:/Users/zoo/Dropbox/SLM/LUTs/SLM_3331_encrypt.txt"; // slm lut dir; change the dir here

unsigned int       n_boards_found = 0U;
bool               constructed_okay = true;



Blink_SDK sdk(bits_per_pixel, pixel_dimension, &n_boards_found,
	&constructed_okay, is_nematic_type, RAM_write_enable,
	use_GPU_if_available, 20U, regional_lut_file);

int getTransformMatrix()
{
	QMessageBox Msgbox;
	MATFile *pmat = matOpen(transform_file, "r");
	if (pmat == NULL) {
		Msgbox.setText("Error reading transform file");
		Msgbox.exec();
		return(1);
	}

	mxArray *parray;
    parray = matGetVariable(pmat, transform_name);
	if (parray == NULL) {
		Msgbox.setText("Error reading transform matrix");
		Msgbox.exec();
		return(1);
	}

	double *pdata;
	pdata = mxGetPr(parray);
	// const size_t *dim_parray = mxGetDimensions(parray);
	for (int i = 0; i < 3; i++)
	{
		Ta[i] = pdata[i];
		Tb[i] = pdata[i + 3];
		Tc[i] = pdata[i + 6];
	}

	matClose(pmat);
	return 0;
}


int setI(char* text)
{
	return 0;
}


float* normalise(float*array, int num_spots)
{
	float sum = 0;
	for (int i = 0; i < num_spots; i++){
		sum += array[i];
	}
	for (int i = 0; i < num_spots; i++){
		array[i] = array[i] / sum;
	}
	return array;
}



std::string process(std::string const& s)
{
	std::string::size_type pos = s.find('.');
	if (pos != std::string::npos)
	{
		return s.substr(0, pos);
	}
	else
	{
		return s;
	}
}

void HoloBlink::weightIntensity()
{
	double distance;
	for (int i = 0; i < par.num_spots; i++)
	{

		distance = sqrt(par.xx[i] * par.xx[i] + par.yy[i] * par.yy[i]);
		par.II[i] = p_weight[1] / (p_weight[0] * SteepnessFudgeFactor*distance + p_weight[1]);

	}
}


int HoloBlink::invertInputs(float*inputs, int num_spots)
{
	// lower input values will get higher weights
//	normalise(inputs, num_spots);
	for (int i = 0; i < num_spots; i++){
		par.I[i] = par.I[i]/inputs[i];

	}
	normalise(par.I, num_spots);
	return 0;
}

float HoloBlink::getTotalPowerOnSample() // not used
{
	float *max_I = std::max_element(par.II.data(), par.II.data() + par.num_spots);
	float totalPower = par.PowerPerCell / *max_I;
	QString str = "Power on sample = " + QString::number(totalPower);
	ui->status_label->setText(str);
	return totalPower;

}

void HoloBlink::acceptConnection()
{
	// need to grab the socket
	client = server->nextPendingConnection();

	client->write("Hello client\r\n");
	client->flush();
	//client->waitForBytesWritten(3000); // not recommended in windows, commented out 20180423

	connect(client, SIGNAL(readyRead()),
		this, SLOT(startRead()));
	ui->status_label->setText("New connection received");
	pLogger->info("New connection received ");
	FLAG_CONNECTED = TRUE;
}

void HoloBlink::startRead()
{
	// keyword and message length -- this need to be the same as RTAOI!
	//ui->status_label->setText("started reading");
	pLogger->info("Started reading");
	if(FLAG_CONNECTED)
	{
	char keyword;
	char * length = new char[4];
	try
	{
		client->read(&keyword, 1);
		client->read(length, 4);
	}
	catch (std::exception& e)
	{
		ui->status_label->setText(e.what());
		pLogger->error(e.what());
		FLAG_CONNECTED = false;
	}

	if (readMessage(keyword, atoi(length)))
	{
		ui->status_label->setText("error reading message");
		pLogger->error("error reading message");
	}

	}

}

int HoloBlink::readMessage(char keyword, int length)
{
	std::clock_t timer_start;
	double duration;
	std::string timer_output;

	timer_start = std::clock();
	char *buffer = new char[length];
	float *values = new float[length];
	client->read(buffer, length);

	char *delim = ";"; // input separated by spaces
	char *pch;
	pch = strtok(buffer, delim);
	int count = 0;
	while (pch != NULL) // need to be less than a byte!

	{
		char *unconverted;
		values[count] = strtof(pch, &unconverted);
		pch = strtok(NULL, delim);
		count += 1;
	}
	duration = (std::clock() - timer_start) / (double)CLOCKS_PER_SEC;
	std::cout << keyword << std::endl;

	timer_output = "recv msg time:" + std::to_string(duration);
	std::cout << timer_output << std::endl;



	// update par
	try{
		int buffer_len = 0;
		timer_start = std::clock();
		switch (keyword)
		{
		case('D') : // reconfigure photostim duration
			par.photoDuration = values[0];
			par.photoSamples = int(par.photoDuration*0.001*par.daqSampleRate) + 1;

			// resize daq buffer
			buffer_len = par.photoSamples * 2;
			WriteData.resize(buffer_len);
			for (int i = 0; i < 5; i++) // 5 ms trigger to spiral
				WriteData[i] = 5;
			WriteData[5] = 0;

			for (int i = par.photoSamples; i < par.photoSamples + 5; i++) // add an offset to avoid burning single point
				WriteData[i] = 0;
			for (int i = par.photoSamples + 5; i < buffer_len - 1; i++) // volt to aom
				WriteData[i] = par.photoVolt;
			WriteData[buffer_len - 1] = 0; // reset to zero

			FLAG_DURATION_UPDATED = TRUE;
			break;

		case('T') : // photostim trigger alone (given aom voltage)
			par.photoVolt = values[0];
			par.photoSamples = int(par.photoDuration*0.001*par.daqSampleRate) + 1;

			// resize daq buffer
			buffer_len = par.photoSamples * 2;
			for (int i = par.photoSamples; i < par.photoSamples+5; i++) // add an offset to avoid burning single point
				WriteData[i] = 0;
			for (int i = par.photoSamples+5; i < buffer_len - 1; i++) // volt to aom
				WriteData[i] = par.photoVolt;
			WriteData[buffer_len - 1] = 0; // reset to zero
			FLAG_TRIGGER_PHOTO = TRUE;
			break;

		case('P') : // coordinates and aom voltage 
			par.num_spots = int((count - 1) / 2);
			par.xx_i.resize(par.num_spots);
			par.yy_i.resize(par.num_spots);
			par.photoVolt = float64(values[count - 1]);

			for (int i = 0; i < par.num_spots; i++)
			{
				par.xx_i[i] = values[i];
				par.yy_i[i] = values[i + par.num_spots];
			}


			// resize daq buffer
			buffer_len = par.photoSamples * 2;
			WriteData.resize(buffer_len);
			for (int i = 0; i < 5; i++) // 5 ms trigger to spiral
				WriteData[i] = 5;
			WriteData[5] = 0;

			for (int i = par.photoSamples; i < par.photoSamples + 5; i++) // add an offset to avoid burning single point
				WriteData[i] = 0;
			for (int i = par.photoSamples + 5; i < buffer_len - 1; i++) // volt to aom
				WriteData[i] = par.photoVolt;
			WriteData[buffer_len - 1] = 0; // reset to zero

			OutputDebugStringA("got new coords and voltage ");
			FLAG_XY_CHANGE = TRUE;
			FLAG_TRIGGER_PHOTO = TRUE;
			break;

		case('C') : // coordinates alone
			par.num_spots = int(count / 2);
			par.xx_i.resize(par.num_spots);
			par.yy_i.resize(par.num_spots);

			for (int i = 0; i < par.num_spots; i++)
			{
				par.xx_i[i] = values[i];
				par.yy_i[i] = values[i + par.num_spots];
			}
			OutputDebugStringA("got new coords");
			FLAG_XY_CHANGE = TRUE;
			break;

		case('I') :
			par.I = values;
			FLAG_ICHANGE = TRUE;
			break;

		default:
			ui->status_label->setText("no change applied");
		}
	}
	catch (std::exception& e)
	{
		
		pLogger->error("parse parameter exception");
		pLogger->error(e.what());
		ui->status_label->setText(e.what());
		return 1;
	}
	duration = (std::clock() - timer_start) / (double)CLOCKS_PER_SEC;
	timer_output = "update parameters time:" + std::to_string(duration);
	std::cout << timer_output << std::endl;


	try{
		if (FLAG_XY_CHANGE && (!FLAG_ICHANGE))
		{
			// uniform weight
			//par.I = new float[par.num_spots];
			//par.z = new float[par.num_spots];

			// resize vectors
			par.II.resize(par.num_spots);
			par.zz.resize(par.num_spots);

			for (int i = 0; i < par.num_spots; i++)
			{
				par.II[i] = 1;
				par.zz[i] = 0;
			}
			transformXY_vect();
			// get new hologram and gui
			timer_start = std::clock();
			get_hologram();
			std::cout << "got hologram" << std::endl;
			duration = (std::clock() - timer_start) / (double)CLOCKS_PER_SEC;
			timer_output = "hologram computation time:" + std::to_string(duration);
			std::cout << timer_output << std::endl;

			timer_start = std::clock();
			on_putOnSLM_pushButton_clicked();
			pLogger->info("hologram updated");

			// write to analogue output
			if (FLAG_TRIGGER_PHOTO)
			{
				int ERROR_FLAG = DAQmxWriteAnalogF64(WriteTaskHandle, par.photoSamples, 0, 10, DAQmx_Val_GroupByChannel, WriteData.data(), &written, NULL);
				Sleep(par.photoDuration);
				FLAG_TRIGGER_PHOTO = FALSE;
				std::cout << "triggers sent" << std::endl;
				pLogger->info("photostim sent");
			}

			client->write("Done\r\n",7); //sizeof("Done\r\n")
			client->flush();
			// making sure echo is written - waitForBytesWritten will slow things down, try if still not working
			// if(!client->waitForBytesWritten(1000)){
			// 	pLogger->info("wait for bytes written timed out")
			// }

			pLogger->info("Returned to client");
			std::cout << "put hologram on SLM, Done" << std::endl;
			duration = (std::clock() - timer_start) / (double)CLOCKS_PER_SEC;
			timer_output = "update SLM display:" + std::to_string(duration);
			std::cout << timer_output << std::endl;

			updateSpots();
			std::cout << "updated spots window " << std::endl;
			updateImages();
			std::cout << "updated hologram window" << std::endl;
			pLogger->info("updated gui");

			FLAG_XY_CHANGE = FALSE;

		}
		if (FLAG_TRIGGER_PHOTO)
		{
			DAQmxWriteAnalogF64(WriteTaskHandle, par.photoSamples, 0, 10, DAQmx_Val_GroupByChannel, WriteData.data(), &written, NULL);
			Sleep(par.photoDuration);
			client->write("Done\r\n");
			std::cout << "triggers sent" << std::endl;
			pLogger->info("photostim sent");
			FLAG_TRIGGER_PHOTO = FALSE;
		}
		if (FLAG_ICHANGE)
		{
			// get new hologram and gui
			get_hologram();
			on_putOnSLM_pushButton_clicked();
			client->write("Done\r\n");
			updateSpots();
			updateImages();
			FLAG_ICHANGE = FALSE;
		}

		if (FLAG_DURATION_UPDATED)
		{
			OutputDebugStringA("photostim duration updated");
			std::cout << par.photoSamples << std::endl;
			client->write("Done\r\n");
			FLAG_DURATION_UPDATED = FALSE;
		}
	}
	catch (std::exception &e)
	{
		
		pLogger->error("holograme update or photostim trigger exception:");
		pLogger->error(e.what());
		ui->status_label->setText(e.what());
		return 1;
	}
	client->flush();
	return 0;


}

HoloBlink::HoloBlink(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::HoloBlink)

{
    ui->setupUi(this);
	m_TableHeader << "x" << "y" << "I";
	ui->ROI_tableWidget->setRowCount(3);
	ui->ROI_tableWidget->setVerticalHeaderLabels(m_TableHeader);
	ui->ROI_tableWidget-> horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
	BLANKIMG.load("C:\\Data\\Zoe\\VSprojects\\HoloBlink\\sources\\BLANKIMG.tif"); // change this to relative path
	server = new QTcpServer(this);

	// whenever a user connects, it will emit signal
	connect(server, SIGNAL(newConnection()),
		this, SLOT(acceptConnection()));

	if (!server->listen(QHostAddress::Any, 8888))
	{
		ui->status_label->setText("Server could not start");
	}
	else
	{
		ui->status_label->setText("Server started");
	}

	// read transform file
	int error = getTransformMatrix();

	// write to two channels 
	DAQmxCreateTask("", &WriteTaskHandle);
	DAQmxCreateAOVoltageChan(WriteTaskHandle, "Dev2/ao5:6", "", 0, 10.0, DAQmx_Val_Volts, NULL);  
	DAQmxCfgSampClkTiming(WriteTaskHandle, "", par.daqSampleRate, DAQmx_Val_Rising, DAQmx_Val_ContSamps, 1000);
	DAQmxStartTask(WriteTaskHandle);

	// resize daq buffer
	int buffer_len = par.photoSamples * 2;
	WriteData.resize(buffer_len);
	for (int i = 0; i < 5; i++) // 5 ms trigger to spiral
		WriteData[i] = 5;
	WriteData[5] = 0;

	for (int i = par.photoSamples; i < buffer_len - 1; i++) // set to zero
		WriteData[i] = 0;
	WriteData[buffer_len - 1] = 0;

	// logger
	pLogger = NULL;
	pLogger = Logger::getInstance();
	pLogger->info("Holoblink logger message");
}
HoloBlink::~HoloBlink()
{
    delete ui;
	server->close();
	if (FLAG_CONNECTED){
		client->close();
	}

	sdk.SLM_power(false);

}
void::HoloBlink::on_powerPerCell_valueChanged()
{
	par.PowerPerCell = ui->powerPerCell_spinBox->value();
}
void::HoloBlink::updateIinTable(){
	ui->ROI_tableWidget->setColumnCount(par.num_spots);
	ui->status_label->setText("number of spots =" + QString::number(par.num_spots));

	for (int i = 0; i < par.num_spots; i++){
		// update table
		ui->ROI_tableWidget->setItem(0, i, new QTableWidgetItem(QString::number(par.xx_i[i])));
		ui->ROI_tableWidget->setItem(1, i, new QTableWidgetItem(QString::number(par.yy_i[i])));
		ui->ROI_tableWidget->setItem(2, i, new QTableWidgetItem(QString::number(par.II[i])));
	}
	// getTotalPowerOnSample();
}
void HoloBlink::on_InitialiseSLM_pushButton_clicked()
{
	// Check that everything started up successfully.
	bool okay = constructed_okay && sdk.Is_slm_transient_constructed();

	if (okay)
	{
		enum { e_n_true_frames = 5 };
		sdk.Set_true_frames(e_n_true_frames);
		sdk.SLM_power(true);
		okay = sdk.Load_linear_LUT(board_number);
	}
	else
	{
		ui->status_label->setText(sdk.Get_last_error_message());

	}

	if (okay)
	{
		ui->status_label->setText("SLM connected");
	}
	else
	{
		ui->status_label->setText("SLM not connected");

	}

}
void HoloBlink::on_putOnSLM_pushButton_clicked()
{
	bool okay = constructed_okay && sdk.Write_overdrive_image(1, par.h_pSLM);  // 12 - 13 ms; clear image on CMOS
	if (okay){
		ui->status_label->setText("Update phase mask succeed");
	}
	else{
		ui->status_label->setText(sdk.Get_last_error_message());
		//puts(sdk.Get_last_error_message());
	}

}
void HoloBlink::on_invert_radioButton_clicked()
{
	par.FlagInvert = !par.FlagInvert;
}
void HoloBlink::on_actionOpen_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QString(),
            tr("Image Files (*.tif);;C++ Files (*.cpp *.h)"));

    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly)) {
            QMessageBox::critical(this, tr("Error"), tr("Could not open file"));
            return;
        }
	// display image
	 QPixmap  img(fileName);
	 QGraphicsScene* scene = new QGraphicsScene();
	 scene->addPixmap(img);

	 ui->textEdit->setText(fileName);
	 file.close();

	 // load image as cv mat
	 Mat image;
	 image = imread(fileName.toStdString());   // Read the file

	 // convert rgb to gray, if necessary
	 Mat gray;

	 if (image.channels() == 3)
	 {
		 cvtColor(image, gray, CV_BGR2GRAY);
	 }
	 else
	 {
		 gray = image;
	 }

	 // binarise
	 Mat bw;
	 threshold(gray, bw, 40, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
	 cv::Mat locations;
	 cv::findNonZero(bw, locations);

	 // get cooridnates
	 const int num_spots = locations.total();
	 x_spots.resize(num_spots);
	 y_spots.resize(num_spots);
	 I.resize(num_spots);
	 ui->ROI_tableWidget->setColumnCount(num_spots);

	 for (int i = 0; i < num_spots; i++) {
		 x_spots[i] = locations.at<Point>(i).x;
		 y_spots[i] = locations.at<Point>(i).y;
		 I[i] = gray.at<uchar>((int)y_spots[i], (int)x_spots[i]);

		 // number ROIs
		 QGraphicsTextItem *text  = scene->addText(QString::number(i + 1));
		 text->setPos(x_spots[i], y_spots[i]);
		 text->setDefaultTextColor(QColor(255,0,0));
	 }


	 // update parameters
	 par.file_name = fileName.toStdString();

	 par.num_spots = num_spots;
	 // do  transform
	 par.x = x_spots.data();
	 par.y = y_spots.data();
	// transformXY(par.x, par.y, Ta, Tb, num_spots);

	 par.xx_i.resize(par.num_spots);
	 par.yy_i.resize(par.num_spots);
	 par.zz.resize(par.num_spots);
	 par.II.resize(par.num_spots);


	 for (int i = 0; i < par.num_spots; i++)
	 {
		 par.xx_i[i] = par.x[i];
		 par.yy_i[i] = par.y[i];
		 par.II[i] = 1;
		 par.zz[i] = 0;
	 }
	 OutputDebugStringA("got new coords");
	 transformXY_vect();

	 // not doing transform (for debug)
	 //for (int i = 0; i < num_spots; i++)
	 //{
		// par.x[i] -= 256;
		// par.y[i] -= 256;
	 //}

	 par.I = I.data();
	 normalise(par.I,num_spots);
	 //
	 get_hologram();

	 // show image
	 ui->Image_graphicsView->setScene(scene);
	 ui->Image_graphicsView->show();

	 // show transformed spots
	 for (int i = 0; i < num_spots; i++){
		 // number ROIs
		 QGraphicsTextItem *text = scene->addText(QString::number(i + 1));
		 text->setPos(par.xx[i]+256, par.yy[i]+256);
		 text->setDefaultTextColor(QColor(0, 255, 0)); }
		 updateIinTable();

    }
}
void HoloBlink::updateSpots()
{
	QGraphicsScene* scene = new QGraphicsScene();
	scene->addPixmap(BLANKIMG);
	for (int i = 0; i < par.num_spots; i++){
		// number ROIs
		QGraphicsTextItem *text = scene->addText(QString::number(i + 1));
		text->setPos(par.xx[i] + 256, par.yy[i] + 256);
		text->setDefaultTextColor(QColor(0, 255, 0));
	}
	// show image
	ui->Image_graphicsView->setScene(scene);
	ui->Image_graphicsView->show();
	// updateIinTable(); // dont bother with table display
}
Mat HoloBlink::make_save_img(float* x, float* y, bool IfNormalise)
{
	Mat image(512, 512, CV_8UC1,Scalar(0));
	if (IfNormalise){
		for (int i = 0; i < par.num_spots; i++){
			image.at<uchar>((int)y[i] + 255, (int)x[i] + 255) = 255;
		}
	}
	else{
		for (int i = 0; i < par.num_spots; i++){
			image.at<uchar>((int)y[i], (int)x[i]) = 255;
		}
	}

	return image;

}
void HoloBlink::on_actionSave_triggered()
{
	// save target image; transfrome image and hologram
	// std::string save_name = process(par.file_name);
	QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
		"/home",
		QFileDialog::ShowDirsOnly
		| QFileDialog::DontResolveSymlinks);

	std::string save_dir = dir.toStdString();
	cv::imwrite(save_dir+"/phasemask.bmp", cv::Mat(512, 512, CV_8U, par.h_pSLM));
	cv::imwrite(save_dir + "/transformed_targets.bmp", make_save_img(par.xx.data(),par.yy.data(),TRUE));
	cv::imwrite(save_dir + "/targets.bmp", make_save_img(par.xx_i.data(), par.yy_i.data(), FALSE));
	ui->status_label->setText(QString::fromStdString("Saved to " + save_dir));
}

void HoloBlink::transformXY(float*x, float*y, double *A, double *B, int len)
{
	// save to input coords
	par.x_i = new float(par.num_spots);
	par.y_i = new float(par.num_spots);

	for (int i = 0; i < len; i++)
	{
		float this_x = x[i] + XOFFSET;
		float this_y = y[i] + YOFFSET;
		
		par.x_i[i] = this_x - 255 ;
		par.y_i[i] = this_y - 255 ;

		double this_z = Tc[0] * this_x + Tc[1] * this_y + Tc[2];

		// 3 pixels off in rtaoi
		x[i] = roundf(A[0] * this_x/this_z + A[1] * this_y/this_z + A[2]/this_z) - 256;
		y[i] = roundf(B[0] * this_x/this_z + B[1] * this_y/this_z + B[2] / this_z) - 256;
		
		
	}
}

void HoloBlink::transformXY_vect()
{
	// resize 
	par.xx.resize(par.num_spots);
	par.yy.resize(par.num_spots);

	//
	for (int i = 0; i < par.num_spots; i++)
	{
		float this_x = par.xx_i[i] + XOFFSET;
		float this_y = par.yy_i[i] + YOFFSET;

		double this_z = Tc[0] * this_x + Tc[1] * this_y + Tc[2];

		// 3 pixels off in rtaoi
		par.xx[i] = roundf(Ta[0] * this_x / this_z + Ta[1] * this_y / this_z + Ta[2] / this_z) - 256;
		par.yy[i] = roundf(Tb[0] * this_x / this_z + Tb[1] * this_y / this_z + Tb[2] / this_z) - 256;


	}


}

int HoloBlink::get_hologram()
{
	// weighting 
	weightIntensity();

//	 int cuda_error = GenerateHologram(par.h_test, par.h_pSLM, par.y, par.x, par.z, par.I, par.num_spots, par.N_iterations, par.h_obtainedAmps, par.method);
	// using vectors
	int cuda_error = GenerateHologram(par.h_test, par.h_pSLM, par.yy.data(), par.xx.data(), par.zz.data(), par.II.data(), par.num_spots, par.N_iterations, par.h_obtainedAmps, par.method);

	return cuda_error;
}

void HoloBlink::updateImages()
{
	cv::Mat cv_holo8 = cv::Mat(512, 512, CV_8U, par.h_pSLM);
	//cv_holo.convertTo(cv_holo8, CV_8U, 0.00390625);

	// display hologram
	QImage img(cv_holo8.data, 512, 512, QImage::Format_Grayscale8);
	img.bits();
	QPixmap imgPM = QPixmap::fromImage(img);
	QGraphicsScene* scene = new QGraphicsScene(); //move it somewhere else
	scene->addPixmap(imgPM);
	ui->Hologram_graphicsView->setScene(scene);
	ui->Hologram_graphicsView->show();

	// save hologram to par - maybe this is what caused the heap error!!
	//std::memcpy(par.h_pSLM, cv_holo8.data, 512 * 512 * sizeof(uchar));
}
void HoloBlink::on_transform_pushButton_clicked()
{

	int cuda_error = GenerateHologram(par.h_test, par.h_pSLM, par.y, par.x, par.z, par.I, par.num_spots, par.N_iterations, par.h_obtainedAmps, par.method);

	updateImages();
}


int::HoloBlink::on_insert_pushButton_clicked(){
	QString text = ui->ROI_textEdit->toPlainText();
	QStringList list = text.split(",");

	int length = list.length();
	if (length != par.num_spots) 
	{
		updateIinTable();
		return 1; 
	}

	vector<float> temp;
	temp.resize(par.num_spots);
	for (int i = 0; i < length; i++){
		temp[i] = list[i].toFloat();
	}

	if (!par.FlagInvert)
	{	
		std::copy(temp.data(), temp.data() + length, par.I);
		normalise(par.I, par.num_spots);
	}
	else
	{
		invertInputs(temp.data(), length);
	}
	
	updateIinTable();
	return 0;
}

/********************************************************************************
** Form generated from reading UI file 'holoblink.ui'
**
** Created by: Qt User Interface Compiler version 5.6.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_HOLOBLINK_H
#define UI_HOLOBLINK_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_HoloBlink
{
public:
    QAction *actionOpen;
    QAction *actionSave;
    QWidget *centralWidget;
    QPushButton *transform_pushButton;
    QGraphicsView *Image_graphicsView;
    QGraphicsView *Hologram_graphicsView;
    QTextEdit *textEdit;
    QTableWidget *ROI_tableWidget;
    QTextEdit *ROI_textEdit;
    QPushButton *Insert_pushButton;
    QRadioButton *invertInput_radioButton;
    QLabel *label;
    QLabel *label_2;
    QLabel *status_label;
    QSpinBox *powerPerCell_spinBox;
    QLabel *label_3;
    QPushButton *putOnSLM_pushButton;
    QPushButton *InitialiseSLM_pushButton;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *HoloBlink)
    {
        if (HoloBlink->objectName().isEmpty())
            HoloBlink->setObjectName(QStringLiteral("HoloBlink"));
        HoloBlink->resize(1139, 881);
        actionOpen = new QAction(HoloBlink);
        actionOpen->setObjectName(QStringLiteral("actionOpen"));
        actionSave = new QAction(HoloBlink);
        actionSave->setObjectName(QStringLiteral("actionSave"));
        centralWidget = new QWidget(HoloBlink);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        transform_pushButton = new QPushButton(centralWidget);
        transform_pushButton->setObjectName(QStringLiteral("transform_pushButton"));
        transform_pushButton->setGeometry(QRect(710, 560, 111, 31));
        Image_graphicsView = new QGraphicsView(centralWidget);
        Image_graphicsView->setObjectName(QStringLiteral("Image_graphicsView"));
        Image_graphicsView->setGeometry(QRect(20, 30, 512, 512));
        Hologram_graphicsView = new QGraphicsView(centralWidget);
        Hologram_graphicsView->setObjectName(QStringLiteral("Hologram_graphicsView"));
        Hologram_graphicsView->setGeometry(QRect(560, 30, 512, 512));
        textEdit = new QTextEdit(centralWidget);
        textEdit->setObjectName(QStringLiteral("textEdit"));
        textEdit->setGeometry(QRect(70, 560, 461, 31));
        ROI_tableWidget = new QTableWidget(centralWidget);
        if (ROI_tableWidget->columnCount() < 10)
            ROI_tableWidget->setColumnCount(10);
        if (ROI_tableWidget->rowCount() < 3)
            ROI_tableWidget->setRowCount(3);
        ROI_tableWidget->setObjectName(QStringLiteral("ROI_tableWidget"));
        ROI_tableWidget->setGeometry(QRect(20, 610, 1051, 111));
        ROI_tableWidget->setRowCount(3);
        ROI_tableWidget->setColumnCount(10);
        ROI_textEdit = new QTextEdit(centralWidget);
        ROI_textEdit->setObjectName(QStringLiteral("ROI_textEdit"));
        ROI_textEdit->setGeometry(QRect(70, 740, 681, 31));
        Insert_pushButton = new QPushButton(centralWidget);
        Insert_pushButton->setObjectName(QStringLiteral("Insert_pushButton"));
        Insert_pushButton->setGeometry(QRect(990, 740, 81, 31));
        invertInput_radioButton = new QRadioButton(centralWidget);
        invertInput_radioButton->setObjectName(QStringLiteral("invertInput_radioButton"));
        invertInput_radioButton->setGeometry(QRect(790, 750, 171, 17));
        label = new QLabel(centralWidget);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(30, 570, 47, 13));
        label_2 = new QLabel(centralWidget);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(30, 750, 47, 13));
        status_label = new QLabel(centralWidget);
        status_label->setObjectName(QStringLiteral("status_label"));
        status_label->setGeometry(QRect(20, 810, 881, 16));
        powerPerCell_spinBox = new QSpinBox(centralWidget);
        powerPerCell_spinBox->setObjectName(QStringLiteral("powerPerCell_spinBox"));
        powerPerCell_spinBox->setGeometry(QRect(650, 561, 42, 31));
        powerPerCell_spinBox->setValue(6);
        label_3 = new QLabel(centralWidget);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(570, 570, 71, 16));
        putOnSLM_pushButton = new QPushButton(centralWidget);
        putOnSLM_pushButton->setObjectName(QStringLiteral("putOnSLM_pushButton"));
        putOnSLM_pushButton->setGeometry(QRect(950, 560, 111, 31));
        InitialiseSLM_pushButton = new QPushButton(centralWidget);
        InitialiseSLM_pushButton->setObjectName(QStringLiteral("InitialiseSLM_pushButton"));
        InitialiseSLM_pushButton->setGeometry(QRect(830, 560, 111, 31));
        HoloBlink->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(HoloBlink);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1139, 21));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        HoloBlink->setMenuBar(menuBar);
        mainToolBar = new QToolBar(HoloBlink);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        HoloBlink->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(HoloBlink);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        HoloBlink->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addAction(actionSave);

        retranslateUi(HoloBlink);
        QObject::connect(transform_pushButton, SIGNAL(clicked()), HoloBlink, SLOT(on_transform_pushButton_clicked()));
        QObject::connect(Insert_pushButton, SIGNAL(clicked()), HoloBlink, SLOT(on_insert_pushButton_clicked()));
        QObject::connect(invertInput_radioButton, SIGNAL(clicked()), HoloBlink, SLOT(on_invert_radioButton_clicked()));
        QObject::connect(powerPerCell_spinBox, SIGNAL(valueChanged(int)), HoloBlink, SLOT(on_powerPerCell_valueChanged()));
        QObject::connect(putOnSLM_pushButton, SIGNAL(clicked()), HoloBlink, SLOT(on_putOnSLM_pushButton_clicked()));
        QObject::connect(InitialiseSLM_pushButton, SIGNAL(clicked()), HoloBlink, SLOT(on_InitialiseSLM_pushButton_clicked()));

        QMetaObject::connectSlotsByName(HoloBlink);
    } // setupUi

    void retranslateUi(QMainWindow *HoloBlink)
    {
        HoloBlink->setWindowTitle(QApplication::translate("HoloBlink", "HoloBlink", 0));
        actionOpen->setText(QApplication::translate("HoloBlink", "Open", 0));
        actionSave->setText(QApplication::translate("HoloBlink", "Save", 0));
        transform_pushButton->setText(QApplication::translate("HoloBlink", "Get PhaseMask", 0));
        Insert_pushButton->setText(QApplication::translate("HoloBlink", "Update ", 0));
        invertInput_radioButton->setText(QApplication::translate("HoloBlink", "Invert (inputs are STA peaks)", 0));
        label->setText(QApplication::translate("HoloBlink", "Path", 0));
        label_2->setText(QApplication::translate("HoloBlink", "Set I", 0));
        status_label->setText(QApplication::translate("HoloBlink", "Status", 0));
        label_3->setText(QApplication::translate("HoloBlink", "Max mW/cell", 0));
        putOnSLM_pushButton->setText(QApplication::translate("HoloBlink", "Put on SLM", 0));
        InitialiseSLM_pushButton->setText(QApplication::translate("HoloBlink", "Initialise SLM", 0));
        menuFile->setTitle(QApplication::translate("HoloBlink", "File", 0));
    } // retranslateUi

};

namespace Ui {
    class HoloBlink: public Ui_HoloBlink {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_HOLOBLINK_H

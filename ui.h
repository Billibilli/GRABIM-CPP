#ifndef UI_H
#define UI_H

#include <QClipboard>
#include <QApplication>
#include <QGridLayout>
#include <QPixmap>
#include <QMainWindow>
#include <QGroupBox>
#include <QRadioButton>
#include <QLabel>
#include <QWidget>
#include <QApplication>
#include <QString>
#include <QPushButton>
#include <QtSvg>
#include <QObject>
#include <QFileDialog>
#include "io.h"
#include "GRABIM.h"
#include <sstream>

class ui :public QMainWindow
{
     Q_OBJECT
public:
    ui();

private:
    QPushButton * SourceFileButton, * LoadFileButton;
    QPushButton * RunButton, * CancelButton;
    QPushButton * GNUplotButton, *QucsSchematicPathButton;
    QPushButton * TopoScriptButton;
    QWidget * centralWidget;
    QSvgWidget *imgWidget;
    QComboBox *Topology, *SearchModeCombo;

    QRadioButton *ArbitraryTopology, *SearchBestTopology;
    QComboBox * LocalOptCombo, *minFUnitsCombo, * maxFUnitsCombo;
    QComboBox * QucsCombo;
    QLabel * LocalOptLabel, *minFLabel, *maxFLabel, *CodeLabel;
    QString  SourceFile, LoadFile;
    QString GNUplot_path, TopoScript_path;
    QLineEdit * minFEdit, * maxFEdit, *ArbitraryTopologyLineEdit;

    QCheckBox * FixedZLCheckbox, * FixedZSCheckbox;
    QLineEdit * FixedZSLineedit, *FixedZLLineedit;

    QLabel *ZLOhmLabel, *ZSOhmLabel;

    QString QucsSchematicPath;

 private slots:

    void ArbitraryTopology_clicked();
    void SearchBestTopology_clicked();
    void go_clicked();
    void cancel_clicked();
    void SourceImpedance_clicked();
    void LoadImpedance_clicked();
    void GNUplotOutput_clicked();
    void FixedZSCheckbox_clicked();
    void FixedZLCheckbox_clicked();
    void TopoScriptButton_clicked();
    void QucsSchPathButton_clicked();

    cx_double getComplexImpedanceFromText(char *);

    double getFreqScale(int);
    QString getTopoScriptPath();
};

#endif // UI_H

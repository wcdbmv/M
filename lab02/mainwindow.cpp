#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent),
	ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	getValues();
}

MainWindow::~MainWindow() {
	delete ui;
}

void MainWindow::getValues() {
	parameters.R   = ui->R_doubleSpinBox->value();
	parameters.Le  = ui->Le_doubleSpinBox->value();
	parameters.Lk  = ui->Lk_doubleSpinBox->value();
	parameters.Ck  = ui->Ck_doubleSpinBox->value();
	parameters.Rk  = ui->Rk_doubleSpinBox->value();
	parameters.Uc0 = ui->Uc0_doubleSpinBox->value();
	parameters.I0  = ui->I0_doubleSpinBox->value();
}

void MainWindow::on_calculatePushButton_clicked() {
	//
}

#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include <QMessageBox>

#include "simulate.hpp"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, ui(new Ui::MainWindow)
{
	ui->setupUi(this);
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::on_simulatePushButton_clicked()
{
	SimulateParams params;

	params.t = ui->tSpinBox->value();
	params.dt = ui->dtSpinBox->value();
	params.n = ui->nSpinBox->value();

	params.op1.t = ui->tOp1SpinBox->value();
	params.op1.dt = ui->dtOp1SpinBox->value();
	params.op2.t = ui->tOp2SpinBox->value();
	params.op2.dt = ui->dtOp2SpinBox->value();
	params.op3.t = ui->tOp3SpinBox->value();
	params.op3.dt = ui->dtOp3SpinBox->value();

	params.comp1.t = ui->tComp1SpinBox->value();
	params.comp2.t = ui->tComp2SpinBox->value();

	const auto result = simulateEvent(params);

	ui->nReturnLineEdit->setText(QString::number(result.n_return));
	ui->pReturnLineEdit->setText(QString::number(result.p_return));
}

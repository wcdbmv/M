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

	params.client.t = ui->tClientSpinBox->value();
	params.client.dt = ui->dtClientSpinBox->value();
	params.client.N = ui->nClientSpinBox->value();

	params.terminal.t = ui->tTerminalSpinBox->value();
	params.terminal.dt = ui->dtTerminalSpinBox->value();
	params.terminal.L = ui->lTerminalSpinBox->value();

	params.window1.t = ui->tWindow1SpinBox->value();
	params.window1.dt = ui->dtWindow1SpinBox->value();
	params.window1.L = ui->lWindow1SpinBox->value();
	params.window1.p_return = ui->pReturnWindow1LineEdit->text().toDouble();

	params.window2.t = ui->tWindow2SpinBox->value();
	params.window2.dt = ui->dtWindow2SpinBox->value();
	params.window2.L = ui->lWindow2SpinBox->value();
	params.window2.p_return = 0;

	params.window3.t = ui->tWindow3SpinBox->value();
	params.window3.dt = ui->dtWindow3SpinBox->value();
	params.window3.L = ui->lWindow3SpinBox->value();
	params.window3.p_return = 0;

	const auto result = simulateEvent(params);

	ui->nDroppedTerminalLineEdit->setText(QString::number(result.terminal.n_dropped));
	ui->pDroppedTerminalLineEdit->setText(QString::number(result.terminal.p_dropped));
	ui->nDroppedWindow1LineEdit->setText(QString::number(result.window1.n_dropped));
	ui->pDroppedWindow1LineEdit->setText(QString::number(result.window1.p_dropped));
	ui->nDroppedWindow2LineEdit->setText(QString::number(result.window2.n_dropped));
	ui->pDroppedWindow2LineEdit->setText(QString::number(result.window2.p_dropped));
	ui->nDroppedWindow3LineEdit->setText(QString::number(result.window3.n_dropped));
	ui->pDroppedWindow3LineEdit->setText(QString::number(result.window3.p_dropped));
}

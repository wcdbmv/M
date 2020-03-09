#include "mainwindow.h"
#include "ui_mainwindow.h"



MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent),
	ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	getValues();

	customPlots = {
		ui->I_customPlot,
		ui->Uc_customPlot,
		ui->Rp_customPlot,
		ui->Ucp_customPlot
	};

	for (auto&& customPlot: customPlots) {
		customPlot->addGraph();
		customPlot->graph(0)->setPen(QPen(Qt::red));
		customPlot->xAxis->setLabel("t, мкс");
		customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
	}

	ui->I_customPlot->yAxis->setLabel("I, А");
	ui->Uc_customPlot->yAxis->setLabel("Uc, В");
	ui->Rp_customPlot->yAxis->setLabel("Rp, Ом");
	ui->Ucp_customPlot->yAxis->setLabel("Ucp, В");
}

MainWindow::~MainWindow() {
	delete ui;
}

void MainWindow::getValues() {
	parameters.R   = ui->R_doubleSpinBox->value();
	parameters.Le  = ui->Le_doubleSpinBox->value();
	parameters.Lk  = ui->Lk_doubleSpinBox->value() * 1e-6;
	parameters.Ck  = ui->Ck_doubleSpinBox->value() * 1e-6;
	parameters.Rk  = ui->Rk_doubleSpinBox->value();
	parameters.Uc0 = ui->Uc0_doubleSpinBox->value();
	parameters.I0  = ui->I0_doubleSpinBox->value();
}

void MainWindow::on_calculatePushButton_clicked() {
	getValues();

	auto dependecy = solve(parameters);

	ui->I_customPlot->graph(0)->setData(dependecy.t, dependecy.I);
	ui->Uc_customPlot->graph(0)->setData(dependecy.t, dependecy.Uc);
	ui->Rp_customPlot->graph(0)->setData(dependecy.t, dependecy.Rp);
	ui->Ucp_customPlot->graph(0)->setData(dependecy.t, dependecy.Ucp);

	for (auto&& customPlot: customPlots) {
		customPlot->graph(0)->rescaleAxes();
		customPlot->replot();
	}
}

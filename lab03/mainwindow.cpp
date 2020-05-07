#include "mainwindow.hpp"
#include "./ui_mainwindow.h"
#include "solve.hpp"

MainWindow::MainWindow(QWidget *parent):
		QMainWindow(parent),
		ui(new Ui::MainWindow) {
	ui->setupUi(this);

	ui->customPlot->addGraph();
	ui->customPlot->graph(0)->setPen(QPen(Qt::red));
	ui->customPlot->xAxis->setLabel("x, см");
	ui->customPlot->yAxis->setLabel("T, К");
	ui->customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

MainWindow::~MainWindow() {
	delete ui;
}

void MainWindow::on_calculatePushButton_clicked() {
	const Parameters parameters = {
		ui->k0_doubleSpinBox->value(),
		ui->kN_doubleSpinBox->value(),
		ui->alpha0_doubleSpinBox->value(),
		ui->alphaN_doubleSpinBox->value(),
		ui->F0_doubleSpinBox->value(),
	};

	const auto dependency = solve(parameters);

	ui->customPlot->graph(0)->setData(dependency.x, dependency.T);
	ui->customPlot->graph(0)->rescaleAxes();
	ui->customPlot->replot();
}

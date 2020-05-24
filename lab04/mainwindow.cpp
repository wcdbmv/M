#include "mainwindow.hpp"
#include "./ui_mainwindow.h"
#include "solve.hpp"

MainWindow::MainWindow(QWidget *parent):
		QMainWindow(parent),
		ui(new Ui::MainWindow) {
	ui->setupUi(this);

	customPlots = {
		ui->Tx_customPlot,
		ui->Tt_customPlot,
	};

	for (auto&& customPlot: customPlots) {
		customPlot->yAxis->setLabel("T, К");
		customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
		customPlot->legend->setVisible(true);
		customPlot->legend->setRowSpacing(-3);
	}

	ui->Tx_customPlot->xAxis->setLabel("x, см");
	ui->Tt_customPlot->xAxis->setLabel("t, мкс");
}

MainWindow::~MainWindow() {
	delete ui;
}

void MainWindow::on_calculatePushButton_clicked() {
	const Parameters parameters = {
		ui->F0_doubleSpinBox->value(),
	};

	auto result = solve(parameters);

	ui->Tx_customPlot->clearGraphs();
	for (int i = 0, j = 0; i < result.Tx.size(); i += 3, ++j) {
		auto label = QString("time = ") + QString::number(result.t[i]);
		ui->Tx_customPlot->addGraph();
		QPen pen(QColor(int(qSin(j*0.3)*100)+100, int(qSin(j*0.6+0.7)*100)+100, int(qSin(j*0.4+0.6)*100)+100));
		pen.setWidth(2);
		ui->Tx_customPlot->graph()->setPen(pen);
		ui->Tx_customPlot->graph()->setName(label);
		ui->Tx_customPlot->graph()->setData(result.x, result.Tx[i]);
	}

	ui->Tt_customPlot->clearGraphs();
	for (int i = 0; i < result.Tt.size(); ++i) {
		auto label = QString("len = ") + QString::number(i);
		ui->Tt_customPlot->addGraph();
		QPen pen(QColor(int(qSin(i*0.3)*100)+100, int(qSin(i*0.6+0.7)*100)+100, int(qSin(i*0.4+0.6)*100)+100));
		pen.setWidth(2);
		ui->Tt_customPlot->graph()->setPen(pen);
		ui->Tt_customPlot->graph()->setName(label);
		ui->Tt_customPlot->graph()->setData(result.t, result.Tt[i]);
	}

	for (auto&& customPlot: customPlots) {
		customPlot->rescaleAxes();
		customPlot->replot();
	}
}

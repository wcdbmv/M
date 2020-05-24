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


using ModelResult = std::pair<QVector<double>, QVector<std::pair<double, QVector<double>>>>;

QVector<double> form_time_array(const ModelResult &results){
	QVector<double> times(results.second.size());
	auto get_time = [](std::pair<double, QVector<double>> item){return item.first;};
	std::transform(results.second.begin(), results.second.end(), times.begin(), get_time);
	return times;
}

QVector<QVector<double>> form_temp_from_time_values(const ModelResult &results){
	QVector<QVector<double>> result(11);
	auto step = static_cast<int>(1. / (results.first[1] - results.first[0]));

	for (auto i = 0; i <= 10; i++) {
		for (auto &item: results.second) {
			result[i].push_back(item.second[i * step]);
		}
	}

	return result;
}

void MainWindow::on_calculatePushButton_clicked() {
	 const Parameters parameters = {
		ui->F0_doubleSpinBox->value(),
	};

	auto result = solve(parameters);

	ui->Tx_customPlot->clearGraphs();
	for (auto i = 0, j = 0; i < result.second.size(); i += 3, ++j) {
		auto label = QString("time = ") + QString::number(result.second[i].first);
		ui->Tx_customPlot->addGraph();
		QPen pen(QColor(int(qSin(j*0.3)*100)+100, int(qSin(j*0.6+0.7)*100)+100, int(qSin(j*0.4+0.6)*100)+100));
		pen.setWidth(2);
		ui->Tx_customPlot->graph()->setPen(pen);
		ui->Tx_customPlot->graph()->setName(label);
		ui->Tx_customPlot->graph()->setData(result.first, result.second[i].second);
	}

	ui->Tt_customPlot->clearGraphs();
	auto times = form_time_array(result);
	auto temp_from_time = form_temp_from_time_values(result);
	for (auto i = 1; i < temp_from_time.size(); i++) {
		auto label = QString("len = ") + QString::number(i);
		ui->Tt_customPlot->addGraph();
		QPen pen(QColor(int(qSin(i*0.3)*100)+100, int(qSin(i*0.6+0.7)*100)+100, int(qSin(i*0.4+0.6)*100)+100));
		pen.setWidth(2);
		ui->Tt_customPlot->graph()->setPen(pen);
		ui->Tt_customPlot->graph()->setName(label);
		ui->Tt_customPlot->graph()->setData(times, temp_from_time[i]);
	}

//	const auto dependency = solve(parameters);

//	ui->Tx_customPlot->graph(0)->setData(dependency.x, dependency.Tx[1]);
//	ui->Tt_customPlot->graph(0)->setData(dependency.t, dependency.Tt[1]);

	for (auto&& customPlot: customPlots) {
		customPlot->rescaleAxes();
		customPlot->replot();
	}
}

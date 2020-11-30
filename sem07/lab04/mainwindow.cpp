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

	params.a = ui->aDoubleSpinBox->value();
	params.b = ui->bDoubleSpinBox->value();
	if (params.b <= params.a) {
		QMessageBox::critical(this, "Ошибка", "`a` должно быть меньше `b`");
		return;
	}

	params.mu = ui->muDoubleSpinBox->value();
	params.sigma = ui->sigmaDoubleSpinBox->value();
	if (params.sigma < 0.0) {
		QMessageBox::critical(this, "Ошибка", "`\\sigma` должна быть больше `0`");
		return;
	}

	params.n = ui->nSpinBox->value();
	params.p_return = ui->pReturnDoubleSpinBox->value();
	params.dt = ui->dtLineEdit->text().toDouble();
	if (params.n < 0) {
		QMessageBox::critical(this, "Ошибка", "`N` должно быть не меньше `0`");
		return;
	}
	if (!(0 <= params.p_return && params.p_return <= 1)) {
		QMessageBox::critical(this, "Ошибка", "`p_{возвр}` должно быть от `0` до `1`");
		return;
	}
	if (params.dt <= 0.0) {
		QMessageBox::critical(this, "Ошибка", "`\\Delta t` должно быть больше `0`");
		return;
	}

	auto resultsDt = simulateDt(params);
	auto resultsEvent = simulateEvent(params);

	ui->nReturnDtLineEdit->setText(QString::number(resultsDt.n_return));
	ui->lMaxDtLineEdit->setText(QString::number(resultsDt.max_queue_size));
	ui->nReturnEventLineEdit->setText(QString::number(resultsEvent.n_return));
	ui->lMaxEventLineEdit->setText(QString::number(resultsEvent.max_queue_size));
}

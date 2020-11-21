#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include <QMessageBox>

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
	const auto a = ui->aDoubleSpinBox->value();
	const auto b = ui->bDoubleSpinBox->value();
	if (b <= a) {
		QMessageBox::critical(this, "Ошибка", "`a` должно быть меньше `b`");
		return;
	}

	const auto mu = ui->muDoubleSpinBox->value();
	const auto sigmaSquared = ui->sigmaSquaredDoubleSpinBox->value();
	if (sigmaSquared < 0.0) {
		QMessageBox::critical(this, "Ошибка", "`\\sigma^2` должна быть больше `0`");
		return;
	}

	const auto N = ui->nSpinBox->value();
	const auto pReturn = ui->pReturnDoubleSpinBox->value();
	const auto dt = ui->dtLineEdit->text().toDouble();
	if (!(0 <= pReturn && pReturn <= 1)) {
		QMessageBox::critical(this, "Ошибка", "`p_{возвр}` должно быть от `0` до `1`");
		return;
	}
	if (dt <= 0.0) {
		QMessageBox::critical(this, "Ошибка", "`\\Delta t` должно быть больше `0`");
		return;
	}
}

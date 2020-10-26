#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include <array>
#include <QMessageBox>
#include "markov.hpp"
#include "random.hpp"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	const auto nStates = ui->nStatesSpinBox->value();

	intensityMatrix = QVector<QVector<double>>(nStates, QVector<double>(nStates));

	for (auto& tableWidget : std::array{ui->intensityMatrixTableWidget, ui->resultTableWidget}) {
		tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
	}

	resizeIntensityMatrixTableWidget(nStates);
	connect(
		ui->nStatesSpinBox,
		static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
		ui->intensityMatrixTableWidget,
		[this](int value) {
			resizeIntensityMatrix(value);
			resizeIntensityMatrixTableWidget(value);
		}
	);

	resizeResultTableWidget(nStates);
	ui->resultTableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);
	ui->resultTableWidget->setHorizontalHeaderLabels(QStringList() << "p" << "t1" << "t2");
}

MainWindow::~MainWindow()
{
	delete ui;
}

#include <QtDebug>
void MainWindow::on_calculatePushButton_clicked()
{
	const auto nStates = ui->nStatesSpinBox->value();

	QVector<double> p0_1(nStates);
	p0_1[0] = 1;
	QVector<double> p0_a(nStates, 1.0 / nStates);

	qDebug() << intensityMatrix;
	const auto result = solve(intensityMatrix);
	const auto time_result_1 = get_system_times(intensityMatrix, result, p0_1);
	const auto time_result_a = get_system_times(intensityMatrix, result, p0_a);

	resizeResultTableWidget(nStates);
	for (int i = 0; i < nStates; ++i) {
		ui->resultTableWidget->item(i, 0)->setText(QString::number(result[i]));
		ui->resultTableWidget->item(i, 1)->setText(QString::number(time_result_1[i]));
		ui->resultTableWidget->item(i, 2)->setText(QString::number(time_result_a[i]));
	}
}

void MainWindow::resizeIntensityMatrix(int size)
{
	intensityMatrix.resize(size);
	for (auto& row : intensityMatrix) {
		row.resize(size);
	}
}

void MainWindow::resizeTableWidget(QTableWidget *tableWidget, int rows, int cols)
{
	tableWidget->setRowCount(rows);
	tableWidget->setColumnCount(cols);

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			if (!tableWidget->item(i, j)) {
				tableWidget->setItem(i, j, new QTableWidgetItem);
			}
		}
	}
}

void MainWindow::resizeIntensityMatrixTableWidget(int size)
{
	resizeTableWidget(ui->intensityMatrixTableWidget, size, size);
}

void MainWindow::resizeResultTableWidget(int size)
{
	resizeTableWidget(ui->resultTableWidget, size, 3);
}

void MainWindow::on_intensityMatrixTableWidget_cellChanged(int row, int column)
{
	auto item = ui->intensityMatrixTableWidget->item(row, column);
	auto& cell = intensityMatrix[row][column];

	bool ok = true;
	const auto text = item->text();
	const auto value =  text != "" ? text.toDouble(&ok) : 0.0;
	if (!ok) {
		QMessageBox::critical(this, "Ошибка", text + " не является числом");
		item->setText(cell != 0.0 ? QString::number(cell) : "");
		return;
	}
	cell = value;
}

void MainWindow::on_generateIntensityMatrixpushButton_clicked()
{
	const auto nStates = ui->nStatesSpinBox->value();
	const auto maxIntensity = ui->maxIntensityDoubleSpinBox->value();
	for (int i = 0; i < nStates; ++i) {
		for (int j = 0; j < nStates; ++j) {
			double value = 0.0;
			if (i != j) {
				value = randreal(0.0, maxIntensity);
			}
			ui->intensityMatrixTableWidget->item(i, j)->setText(QString::number(value));
		}
	}
}

#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include <array>
#include <QMessageBox>
#include <QThread>
#include "markov.hpp"

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
}

MainWindow::~MainWindow()
{
	delete ui;
}

#include <QtDebug>
void MainWindow::on_calculatePushButton_clicked()
{
	const auto nStates = ui->nStatesSpinBox->value();

	qDebug() << intensityMatrix;
	const auto time_result = get_system_times(intensityMatrix);

	ui->resultTableWidget->setRowCount(nStates);
	for (int i = 0; i < nStates; ++i) {
		ui->resultTableWidget->setItem(i, 0, new QTableWidgetItem(QString::number(time_result[i])));
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
	resizeTableWidget(ui->resultTableWidget, size, 1);
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

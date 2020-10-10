#include "mainwindow.hpp"
#include "./ui_mainwindow.h"

#include "random.hpp"
#include "random_number_table.hpp"

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

void MainWindow::on_generateRandomNumbersWithComputationalMethodPushButton_clicked()
{
	const auto count = ui->amountOfNumbersToGenerateWithComputationalMethod->value();
	if (count <= 0) {
		return;
	}

	const auto ones = generate_sequence<QVector<int>>(count, []() { return randint(0, 9); });
	const auto tens = generate_sequence<QVector<int>>(count, []() { return randint(10, 99); });
	const auto hundreds = generate_sequence<QVector<int>>(count, []() { return randint(100, 999); });

	ui->computationalMethodRandomNumbersTableWidget->fill(ones, tens, hundreds);
}

void MainWindow::on_generateRandomNumbersWithTableMethodPushButton_clicked()
{
	const auto count = ui->amountOfNumbersToGenerateWithTableMethod->value();
	if (count <= 0) {
		return;
	}

	static RandomNumberTable randomNumberTable("../lab01/digits.txt");

	const auto ones = generate_sequence<QVector<int>>(count, []() { return limit(randomNumberTable(), 0, 9); });
	const auto tens = generate_sequence<QVector<int>>(count, []() { return limit(randomNumberTable(), 10, 99); });
	const auto hundreds = generate_sequence<QVector<int>>(count, []() { return limit(randomNumberTable(), 100, 999); });

	ui->tableMethodRandomNumbersTableWidget->fill(ones, tens, hundreds);
}

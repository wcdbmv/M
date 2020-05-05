#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <array>
#include <QMainWindow>
#include "qcustomplot.h"
#include "solve.hpp"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
	Q_OBJECT

public:
	explicit MainWindow(QWidget *parent = nullptr);
	~MainWindow();

private slots:
	void on_calculatePushButton_clicked();

private:
	void getValues();

private:
	Ui::MainWindow *ui;
	std::array<QCustomPlot*, 4> customPlots;

	Parameters parameters;
};

#endif // MAINWINDOW_H

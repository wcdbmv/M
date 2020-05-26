#ifndef MAINWINDOW_HPP_
#define MAINWINDOW_HPP_

#include <array>
#include <QMainWindow>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = nullptr);
	~MainWindow();

private slots:
	void on_calculatePushButton_clicked();

private:
	Ui::MainWindow *ui;
	std::array<QCustomPlot*, 2> customPlots;
};

#endif  // MAINWINDOW_HPP_

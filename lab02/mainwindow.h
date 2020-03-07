#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "parameters.hpp"

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

	Parameters parameters;
};

#endif // MAINWINDOW_H

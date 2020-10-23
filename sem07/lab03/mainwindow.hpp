#pragma once

#include <QMainWindow>
#include <QTableWidget>

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

	void on_intensityMatrixTableWidget_cellChanged(int row, int column);

private:
	Ui::MainWindow *ui;

	QVector<QVector<double>> intensityMatrix;

	void resizeIntensityMatrix(int size);
	void resizeTableWidget(QTableWidget *tableWidget, int rows, int cols);
	void resizeIntensityMatrixTableWidget(int size);
	void resizeResultTableWidget(int size);
};

#include "mainwindow.hpp"
#include "model.hpp"
#include "./ui_mainwindow.h"
#include "plot.hpp"

MainWindow::MainWindow(QWidget* parent) :
	QMainWindow(parent),
	ui(new Ui::MainWindow)
{
	ui->setupUi(this);

	Plot(ui->chart_AB).setChart();
	Plot(ui->chart_Impulse).setChart();
	Plot(ui->chart_Xn).setChart();
}

MainWindow::~MainWindow() {
	delete ui;
}

#include "QtDebug"

void MainWindow::on_calculatePushButton_clicked() {
	const Parameters parameters = {
		ui->Fmax_lineEdit->text().toDouble(),
		ui->Tmax_lineEdit->text().toDouble(),
		ui->nu_lineEdit->text().toDouble(),
		ui->table_h_checkBox->isChecked(),
		ui->table_tau_checkBox->isChecked(),
		ui->plot_T_checkBox->isChecked(),
		ui->plot_c_checkBox->isChecked(),
		ui->plot_impulse_checkBox->isChecked(),
	};

	Model math(parameters);

	qDebug() << "done";

	if (parameters.plot_c) {
		Plot builderAB(ui->chart_AB);

		for (auto graph : math.testAB) {
			builderAB.addGraph();
			double t = 0;
			for (auto el : graph) {
				builderAB.addPoint(t, el);
				t += math.tau_;
			}
			builderAB.editLabel("empty");
		}

		builderAB.editAsixLabels("Время, с", "Температура, K");
		builderAB.setChart();
	}

	if (parameters.plot_impulse) {
		Plot builderImpulse(ui->chart_Impulse);
		builderImpulse.addGraph();

		double t = 0;
		for (auto el : math.testImpulse) {
			builderImpulse.addPoint(t, el);
			t += math.tau_;
		}

		builderImpulse.editLabel("empty");
		builderImpulse.editAsixLabels("Время, с", "Температура, K");
		builderImpulse.setChart();
	}

	if (parameters.plot_T) {
		Plot builderXn(ui->chart_Xn);

		double x = 0;
		for (int i = 0; i < math.temp[0].count() - 1; i += 10, x += 10 * math.h_) {
			builderXn.addGraph();
			double t = 0;
			for (int j = 0; j < math.temp.count(); ++j, t += math.tau_) {
				builderXn.addPoint(t, math.temp[j][i]);
			}
			builderXn.editLabel(QString::number(x));
		}

		builderXn.editAsixLabels("Время, с", "Температура, K");
		builderXn.setChart();
	}

	if (parameters.table_h) {
		QStandardItemModel *model = new QStandardItemModel;
		QStandardItem *item;

		QStringList horizontalHeader;
		horizontalHeader.append("1");
		horizontalHeader.append("0.1");
		horizontalHeader.append("0.01");
		horizontalHeader.append("0.001");

		QStringList verticalHeader;
		for (int i = 0; i < math.testH[0].count(); ++i) {
			verticalHeader.append(QString::number(i * 0.1));

			item = new QStandardItem(QString::number(math.testH[0][i]));
			model->setItem(i, 0, item);

			item = new QStandardItem(QString::number(math.testH[1][i]));
			model->setItem(i, 1, item);

			item = new QStandardItem(QString::number(math.testH[2][i]));
			model->setItem(i, 2, item);

			item = new QStandardItem(QString::number(math.testH[3][i]));
			model->setItem(i, 3, item);
		}

		model->setHorizontalHeaderLabels(horizontalHeader);
		model->setVerticalHeaderLabels(verticalHeader);

		ui->tableView_x->setModel(model);

		ui->tableView_x->resizeRowsToContents();
		ui->tableView_x->resizeColumnsToContents();
	}

	if (parameters.table_tau) {
		QStandardItemModel *model = new QStandardItemModel;
		QStandardItem *item;

		QStringList horizontalHeader;
		horizontalHeader = QStringList();
		horizontalHeader.append("1");
		horizontalHeader.append("0.1");
		horizontalHeader.append("0.01");
		horizontalHeader.append("0.001");

		QStringList verticalHeader;
		verticalHeader = QStringList();
		for (int i = 0; i < math.testTau[0].count(); ++i) {
			verticalHeader.append(QString::number(i * 0.01));

			item = new QStandardItem(QString::number(math.testTau[0][i]));
			model->setItem(i, 0, item);

			item = new QStandardItem(QString::number(math.testTau[1][i]));
			model->setItem(i, 1, item);

			item = new QStandardItem(QString::number(math.testTau[2][i]));
			model->setItem(i, 2, item);

			item = new QStandardItem(QString::number(math.testTau[3][i]));
			model->setItem(i, 3, item);
		}

		model->setHorizontalHeaderLabels(horizontalHeader);
		model->setVerticalHeaderLabels(verticalHeader);

		ui->tableView_t->setModel(model);

		ui->tableView_t->resizeRowsToContents();
		ui->tableView_t->resizeColumnsToContents();
	}
}

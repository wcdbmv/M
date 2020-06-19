#ifndef PLOT_HPP_
#define PLOT_HPP_

#include <QString>
#include <QtCharts>
#include <QChart>

using namespace QtCharts;

class Plot {
public:
	Plot(QChartView *view);
	void addGraph();
	void addPoint(double x, double y);
	void editLabel(QString label);
	void editAsixLabels(QString x, QString y);
	void setChart();

private:
	QChart *_chart;
	QChartView *_view;
	QVector<QLineSeries *> _series;
	QChart::ChartTheme chartTheme = QChart::ChartThemeLight;
};

#endif  // PLOT_HPP_

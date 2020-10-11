#pragma once

#include <QTableWidget>
#include <QHeaderView>

class NumbersTableWidget : public QTableWidget
{
	Q_OBJECT
public:
	explicit NumbersTableWidget(QWidget *parent = nullptr)
		: QTableWidget(parent)
	{
		setColumnCount(3);
		setHorizontalHeaderLabels(QStringList() << "0..9" << "10..99" << "100..999");

		horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
		setEditTriggers(QAbstractItemView::NoEditTriggers);
	}

	void clearRows()
	{
		clearContents();
		model()->removeRows(0, rowCount());
	}

	template <typename T>
	void fill(const QVector<T> &ones, const QVector<T> &tens, const QVector<T> &hundreds)
	{
		const auto rowCount = ones.size();
		Q_ASSERT(tens.size() == rowCount && hundreds.size() == rowCount);

		clearRows();
		setRowCount(rowCount);
		Q_ASSERT(columnCount() == 3);

		for (int i = 0; i < rowCount; ++i) {
			setItem(i, 0, new QTableWidgetItem{QString::number(ones[i])});
			setItem(i, 1, new QTableWidgetItem{QString::number(tens[i])});
			setItem(i, 2, new QTableWidgetItem{QString::number(hundreds[i])});
		}
	}
};


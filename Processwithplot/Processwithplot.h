#ifndef PROCESSWITHPLOT_H
#define PROCESSWITHPLOT_H

#include <QWidget>
#include <QVBoxLayout>
#include "qcustomplot.h"

class ProcessWithPlot : public QWidget {
    Q_OBJECT

public:
    explicit ProcessWithPlot(QWidget* parent = nullptr);

private:
    QCustomPlot* customPlot;

    void performFFTAndPlot();
};

#endif // PROCESSWITHPLOT_H

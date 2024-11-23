#include <QApplication>
#include "Processwithplot.h"

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    // Create the widget
    ProcessWithPlot processWithPlot;
    processWithPlot.show();

    return app.exec();
}

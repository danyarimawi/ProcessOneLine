#include "Processwithplot.h"
#include <fftw3.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

ProcessWithPlot::ProcessWithPlot(QWidget *parent) : QWidget(parent) {
    QVBoxLayout *layout = new QVBoxLayout(this);
    customPlot = new QCustomPlot(this);
    layout->addWidget(customPlot);

    // Perform FFT and plot the result
    performFFTAndPlot();
}

void ProcessWithPlot::performFFTAndPlot() {
    // Load data from file
    ifstream file("C:/Users/danya/Downloads/sin_noise_example.csv");
    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        return;
    }

    vector<double> data;
    string line;
    while (getline(file, line)) {
        const auto pos = line.find(',');
        if (pos != string::npos) {
            line.erase(0, pos + 1); // Remove first column
        }
        data.push_back(stod(line)); // Parse second column
    }
    file.close();

    // Compute wavelengths
    const vector<double> c = {645.9716187, 0.15915893, -4.296465E-06, -3.793865E-10};
    const size_t dataSize = data.size();
    vector<double> wavelengths(dataSize);
    for (size_t i = 0; i < dataSize; ++i) {
        wavelengths[i] = c[0] + c[1] * i + c[2] * i * i + c[3] * i * i * i;
    }

    // Interpolation
    double xwavmax = ceil(*max_element(wavelengths.begin(), wavelengths.end()));
    double xwavmin = floor(*min_element(wavelengths.begin(), wavelengths.end()));

    vector<double> new_x(dataSize);
    for (size_t i = 0; i < dataSize; ++i) {
        new_x[i] = static_cast<double>(xwavmin + i * (xwavmax - xwavmin) / (dataSize - 1));
    }

    auto interpolate = [](const vector<double> &x, const vector<double> &y, const vector<double> &new_x) {
        vector<double> interpolated_values(new_x.size());
        for (size_t i = 0; i < new_x.size(); ++i) {
            double x_value = new_x[i];
            auto it = lower_bound(x.begin(), x.end(), x_value);

            if (it == x.end() || it == x.begin()) {
                interpolated_values[i] = 0.0;
                continue;
            }

            size_t index = it - x.begin();
            double x1 = x[index - 1], x2 = x[index];
            double y1 = y[index - 1], y2 = y[index];
            interpolated_values[i] = y1 + (y2 - y1) * (x_value - x1) / (x2 - x1);
        }
        return interpolated_values;
    };

    vector<double> interpIntensities = interpolate(wavelengths, data, new_x);

    // Zero-padding for FFT
    size_t fftSize = 1 << static_cast<int>(ceil(log2(dataSize))); // Power of 2 size
    vector<double> fftInput(fftSize, 0.0);
    for (size_t i = 0; i < interpIntensities.size(); ++i) {
        fftInput[i] = interpIntensities[i];
    }

    // Perform FFT
    fftw_complex *in = fftw_alloc_complex(fftSize);
    fftw_complex *out = fftw_alloc_complex(fftSize);

    for (size_t i = 0; i < fftSize; ++i) {
        in[i][0] = (i < interpIntensities.size()) ? fftInput[i] : 0.0; // Zero-padding
        in[i][1] = 0.0;
    }

    fftw_plan plan = fftw_plan_dft_1d(fftSize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Compute magnitudes and frequencies
    vector<double> magnitudes(fftSize / 2);
    vector<double> frequencies(fftSize / 2);
    double samplingRate = 1.0 / (new_x[1] - new_x[0]); // Sampling rate from wavelengths

    for (size_t i = 0; i < fftSize / 2; ++i) {
        magnitudes[i] = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        frequencies[i] = static_cast<double>(i) * samplingRate / fftSize;
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    // Convert to QVector for plotting
    QVector<double> qFrequencies;
    QVector<double> qMagnitudes;

    for (const auto &freq : frequencies) {
        qFrequencies.append(freq);
    }
    for (const auto &mag : magnitudes) {
        qMagnitudes.append(mag);
    }

    // Plot the FFT results
    customPlot->addGraph();
    customPlot->graph(0)->setData(qFrequencies, qMagnitudes);
    customPlot->xAxis->setLabel("Frequency (Hz)");
    customPlot->yAxis->setLabel("Magnitude");
    customPlot->xAxis->setRange(0, frequencies.back());
    customPlot->yAxis->setRange(0, *max_element(magnitudes.begin(), magnitudes.end()));
    customPlot->replot();
}

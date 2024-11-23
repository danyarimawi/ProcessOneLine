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

    QVector<double> data;
    string line;
    while (getline(file, line)) {
        const auto pos = line.find(',');
        if (pos != string::npos) {
            line.erase(0, pos + 1); // Remove first column
        }
        data.push_back(stod(line)); // Parse second column
    }
    file.close();


    QVector<double> calibratedIntensities = data;

    // Compute wavelengths
    const vector<double> c = { 645.9716187, 0.15915893, -4.296465E-06, -3.793865E-10 };
    const size_t numPixels = data.size();
  

    vector<double> wavelengths(numPixels);
    for (size_t i = 0; i < numPixels; ++i) {
        wavelengths[i] = c[0] + c[1] * i + c[2] * i * i + c[3] * i * i * i;
    }

    // Interpolation
    int interpSize = numPixels * 2; // Double the sampling frequency
    QVector<double> interpIntensities(interpSize);
    QVector<double> interpWavelengths(interpSize);

    for (int i = 0; i < interpSize; ++i) {
        double t = static_cast<double>(i) / (interpSize - 1);
        double origIndex = t * (numPixels - 1);
        int idx = static_cast<int>(origIndex);
        double frac = origIndex - idx;

        if (idx + 1 < numPixels) {
            interpIntensities[i] = (1 - frac) * calibratedIntensities[idx] + frac * calibratedIntensities[idx + 1];
            interpWavelengths[i] = (1 - frac) * wavelengths[idx] + frac * wavelengths[idx + 1];
        }
        else {
            interpIntensities[i] = calibratedIntensities[idx];
            interpWavelengths[i] = wavelengths[idx];
        }
    }

    // Zero-padding for FFT
    size_t fftSize = 1 << static_cast<int>(ceil(log2(interpSize))); // Power of 2 size
    vector<double> fftInput(fftSize, 0.0);
    for (size_t i = 0; i < interpIntensities.size(); ++i) {
        fftInput[i] = interpIntensities[i];
    }

    // Perform FFT
    fftw_complex* in = fftw_alloc_complex(fftSize);
    fftw_complex* out = fftw_alloc_complex(fftSize);

    for (size_t i = 0; i < fftSize; ++i) {
        in[i][0] = fftInput[i]; // Zero-padding
        in[i][1] = 0.0;
    }

    fftw_plan plan = fftw_plan_dft_1d(fftSize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Compute magnitudes and frequencies
    vector<double> magnitudes(fftSize / 2);
    vector<double> frequencies(fftSize / 2);
    double dt = 0.001; // Adjust this based on your actual sampling rate
    double samplingRate = 1.0 / dt; // Sampling rate in Hz


    for (size_t i = 0; i < fftSize / 2; ++i) {
        magnitudes[i] = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        frequencies[i] = static_cast<double>(i) * samplingRate / fftSize;
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

  
    // Plot the FFT results
    QVector<double> qFrequencies(frequencies.size());
    QVector<double> qMagnitudes(magnitudes.size());

    // Copy elements from std::vector to QVector
    std::copy(frequencies.begin(), frequencies.end(), qFrequencies.begin());
    std::copy(magnitudes.begin(), magnitudes.end(), qMagnitudes.begin());

    customPlot->addGraph();
    customPlot->graph(0)->setData(qFrequencies, qMagnitudes);
    customPlot->xAxis->setLabel("Frequency (Hz)");
    customPlot->yAxis->setLabel("Magnitude");
    customPlot->xAxis->setRange(0, frequencies.back());
    customPlot->yAxis->setRange(0, *max_element(magnitudes.begin(), magnitudes.end()));
    customPlot->replot();
}

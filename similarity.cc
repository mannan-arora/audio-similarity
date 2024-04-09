#include "similarity.h"

#include <fftw3.h>
#include <sndfile.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

Similarity::Similarity(std::string originalFilePath, std::string compareFilePath, float sampleRate) {
    this->originalFilePath = originalFilePath;
    this->compareFilePath = compareFilePath;
    this->sampleRate = sampleRate;

    loadAudioFiles();
}

bool Similarity::loadAudioFiles() {
    originalFile.load(originalFilePath);
    compareFile.load(compareFilePath);

    if (originalFile.getSampleRate() != sampleRate) {
        originalFile.setSampleRate(sampleRate);
    }
    if (compareFile.getSampleRate() != sampleRate) {
        compareFile.setSampleRate(sampleRate);
    }

    return true;
}

float Similarity::zcrSimilarity() {
    float originalZCR = computeZCR(originalFile.samples[0]);
    float compareZCR = computeZCR(compareFile.samples[0]);

    return 1 - std::abs(originalZCR - compareZCR);
}

float Similarity::computeZCR(std::vector<float> buffer) {
    int n = buffer.size();
    double zcrSum = 0.0;
    for (int sample = 1; sample < n; ++sample) {
        float currentSample = buffer[sample];
        float previousSample = buffer[sample - 1];

        if ((currentSample >= 0 && previousSample < 0) ||
            (currentSample < 0 && previousSample >= 0)) {
            zcrSum += 1.0;
        }
    }

    float totalSamples = static_cast<float>(n);
    float zcrNormalized = zcrSum / totalSamples;

    return zcrNormalized;
}

float Similarity::rhythmSimilarity() {
    std::vector<float> originalOnsets = energyDifference(originalFile.samples[0]);
    std::vector<float> compareOnsets = energyDifference(compareFile.samples[0]);

    float pearsonCoefficient = pearsonCorrelation(originalOnsets, compareOnsets);

    return 0.5f * (1 + pearsonCoefficient);
}

float Similarity::pearsonCorrelation(std::vector<float> x, std::vector<float> y) {
    if (x.size() != y.size()) {
        return 0.0f;
    }

    float sumX = 0.0f, sumY = 0.0f, sumXY = 0.0f, sumX2 = 0.0f, sumY2 = 0.0f;
    int n = std::min(x.size(), y.size());

    for (int i = 0; i < n; ++i) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
        sumY2 += y[i] * y[i];
    }

    float numerator = sumXY - (sumX * sumY / n);
    float denominator = sqrt((sumX2 - (sumX * sumX / n)) * (sumY2 - (sumY * sumY / n)));

    return numerator / denominator;
}

std::vector<float> Similarity::energyDifference(std::vector<float> buffer) {
    std::vector<float> onsets;

    int frameSize = 1024;
    int hopSize = 512;
    float threshold = 0.5;

    for (size_t i = 0; i < buffer.size() - frameSize; i += hopSize) {
        float energyCurrent = 0.0;
        float energyNext = 0.0;

        for (int j = 0; j < frameSize; ++j) {
            energyCurrent += buffer[i + j] * buffer[i + j];
        }

        for (int j = 0; j < frameSize; ++j) {
            energyNext += buffer[i + j + hopSize] * buffer[i + j + hopSize];
        }

        if ((energyNext - energyCurrent) / energyCurrent > threshold) {
            onsets.push_back(i + hopSize);
        }
    }

    return onsets;
}

std::vector<std::vector<float>> Similarity::computeSpectrogram(const std::vector<float>& audio, int nfft, int hopLength) {
    int numFrames = (audio.size() - nfft) / hopLength + 1;
    std::vector<std::vector<float>> spectrogram(numFrames, std::vector<float>(nfft / 2 + 1, 0.0));

    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (nfft / 2 + 1));
    double* in = (double*)fftw_malloc(sizeof(double) * nfft);
    fftw_plan p = fftw_plan_dft_r2c_1d(nfft, in, out, FFTW_ESTIMATE);

    for (int i = 0; i < numFrames; ++i) {
        for (int j = 0; j < nfft; ++j) {
            in[j] = static_cast<double>(audio[i * hopLength + j]);
        }
        fftw_execute(p);
        for (int j = 0; j < nfft / 2 + 1; ++j) {
            spectrogram[i][j] = static_cast<float>(std::sqrt(out[j][0] * out[j][0] + out[j][1] * out[j][1]));
        }
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return spectrogram;
}

std::vector<float> Similarity::computeSpectralContrast(const std::vector<std::vector<float>>& spectrogram, int numBands) {
    std::vector<float> spectralContrast(spectrogram.size(), 0.0f);

    for (int i = 0; i < spectrogram.size(); ++i) {
        std::vector<double> sortedSpectrum(spectrogram[i].begin(), spectrogram[i].end());
        std::sort(sortedSpectrum.begin(), sortedSpectrum.end());

        int bandSize = sortedSpectrum.size() / numBands;
        for (int j = 0; j < numBands; ++j) {
            double meanValley = std::accumulate(sortedSpectrum.begin() + j * bandSize, sortedSpectrum.begin() + (j + 1) * bandSize, 0.0) / bandSize;
            double meanPeak = std::accumulate(sortedSpectrum.end() - (j + 1) * bandSize, sortedSpectrum.end() - j * bandSize, 0.0) / bandSize;
            spectralContrast[i] += static_cast<float>(meanPeak - meanValley);
        }
        spectralContrast[i] /= numBands;
    }

    return spectralContrast;
}

float Similarity::computeSimilarityScore(const std::vector<float>& contrast1, const std::vector<float>& contrast2) {
    double similarity = 0.0;
    for (int i = 0; i < contrast1.size(); ++i) {
        similarity += std::abs(contrast1[i] - contrast2[i]);
    }
    similarity = 1 - similarity / contrast1.size();
    return similarity;
}

float Similarity::spectralContrastSimilarity() {
    int nfft = 2048;
    int hopLength = 512;
    int numBands = 6;

    auto spectrogram1 = computeSpectrogram(originalFile.samples[0], nfft, hopLength);
    auto spectrogram2 = computeSpectrogram(compareFile.samples[0], nfft, hopLength);

    auto contrast1 = computeSpectralContrast(spectrogram1, numBands);
    auto contrast2 = computeSpectralContrast(spectrogram2, numBands);

    return computeSimilarityScore(contrast1, contrast2);
}

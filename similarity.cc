#include "similarity.h"
#include <cmath>

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

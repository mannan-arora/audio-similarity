#include "similarity.h"

#include <algorithm>
#include <cmath>
#include <numeric>

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

std::vector<std::vector<float>> Similarity::computeSpectrogram(const std::vector<float>& signal, int windowSize, int hopSize) {
    int numFrames = (signal.size() - windowSize) / hopSize + 1;
    std::vector<std::vector<float>> spectrogram(numFrames, std::vector<float>(windowSize / 2 + 1, 0.0f));

    for (int frame = 0; frame < numFrames; ++frame) {
        int start = frame * hopSize;
        for (int bin = 0; bin < windowSize / 2 + 1; ++bin) {
            float sum = 0.0f;
            for (int n = 0; n < windowSize; ++n) {
                float windowValue = 0.54f - 0.46f * cos(2 * M_PI * n / (windowSize - 1));
                sum += windowValue * signal[start + n] * cos(2 * M_PI * bin * n / windowSize);
            }
            spectrogram[frame][bin] = sum * sum;
        }
    }
    return spectrogram;
}

std::vector<float> Similarity::computeSpectralContrast(const std::vector<std::vector<float>>& spectrogram, int numBands) {
    int numFrames = spectrogram.size();
    int freqBinsPerBand = (spectrogram[0].size() - 2) / numBands;
    std::vector<float> contrast(numBands, 0.0f);

    for (int band = 0; band < numBands; ++band) {
        int startBin = band * freqBinsPerBand + 1;
        int endBin = startBin + freqBinsPerBand;
        std::vector<float> bandValues;
        for (int frame = 0; frame < numFrames; ++frame) {
            bandValues.insert(bandValues.end(), spectrogram[frame].begin() + startBin, spectrogram[frame].begin() + endBin);
        }
        std::nth_element(bandValues.begin(), bandValues.begin() + bandValues.size() / 10, bandValues.end());
        float lowerQuantileMean = std::accumulate(bandValues.begin(), bandValues.begin() + bandValues.size() / 10, 0.0f) / (bandValues.size() / 10);
        std::nth_element(bandValues.begin(), bandValues.end() - bandValues.size() / 10, bandValues.end());
        float upperQuantileMean = std::accumulate(bandValues.end() - bandValues.size() / 10, bandValues.end(), 0.0f) / (bandValues.size() / 10);
        contrast[band] = std::log(upperQuantileMean + 1e-6) - std::log(lowerQuantileMean + 1e-6);
    }

    return contrast;
}

float Similarity::spectralContrastSimilarity() {
    auto originalSpectrogram = computeSpectrogram(originalFile.samples[0], 1024, 512);
    auto compareSpectrogram = computeSpectrogram(compareFile.samples[0], 1024, 512);

    auto originalContrast = computeSpectralContrast(originalSpectrogram, 6);
    auto compareContrast = computeSpectralContrast(compareSpectrogram, 6);

    float similarity = 0.0f;
    for (size_t i = 0; i < originalContrast.size(); ++i) {
        similarity += std::pow(originalContrast[i] - compareContrast[i], 2);
    }
    similarity = std::sqrt(similarity);

    return similarity;
}
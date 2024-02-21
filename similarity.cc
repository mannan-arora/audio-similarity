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
    int numSamples = buffer.size();
    double zcrSum = 0.0;
    for (int sample = 1; sample < numSamples; ++sample) {
        float currentSample = buffer[sample];
        float previousSample = buffer[sample - 1];

        if ((currentSample >= 0 && previousSample < 0) ||
            (currentSample < 0 && previousSample >= 0)) {
            zcrSum += 1.0;
        }
    }

    float totalSamples = static_cast<float>(numSamples);
    float zcrNormalized = zcrSum / totalSamples;

    return zcrNormalized;
}
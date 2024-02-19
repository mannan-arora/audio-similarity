#include "similarity.h"

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
    float originalZCR = computeZCR(originalFile.samples);
    float compareZCR = computeZCR(compareFile.samples);

    return 1 - std::abs(originalZCR - compareZCR);
}

float Similarity::computeZCR(std::vector<std::vector<float>> buffer) {
    int numChannels = buffer.size();
    int numSamplesPerChannel = buffer[0].size();
    double zcrSum = 0.0;

    for (int channel = 0; channel < numChannels; ++channel) {
        for (int sample = 1; sample < numSamplesPerChannel; ++sample) {
            float currentSample = buffer[channel][sample];
            float previousSample = buffer[channel][sample - 1];

            // Detect zero crossing by checking the sign change
            if ((currentSample >= 0 && previousSample < 0) ||
                (currentSample < 0 && previousSample >= 0)) {
                zcrSum += 1.0; // Increment the ZCR count
            }
        }
    }

    // Normalize ZCR by dividing by the total number of samples
    double totalSamples = static_cast<double>(numChannels * numSamplesPerChannel);
    double zcrNormalized = zcrSum / totalSamples;

    return zcrNormalized;
}
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

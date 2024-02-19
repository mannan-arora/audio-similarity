#ifndef _SIMILARITY_H_
#define _SIMILARITY_H_

#include "AudioFile.h"
#include <string>

class Similarity {
private:
    std::string originalFilePath;
    std::string compareFilePath;
    float sampleRate;

    AudioFile<float> originalFile;
    AudioFile<float> compareFile;

public:
    bool loadAudioFiles();

    float zcrSimilairty() const;
    float rhythmSimilarity() const;
    float chromaSimilarity() const;
    float energyEnvelopeSimilarity() const;
    float spectralContrastSimilarity() const;
    float perceptualSimilarity() const;

    Similarity(std::string originalFilePath, std::string compareFilePath, float sampleRate = 44100.f);
};

#endif
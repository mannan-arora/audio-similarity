#ifndef _SIMILARITY_H_
#define _SIMILARITY_H_

#include "AudioFile.h"
#include <string>

class Similarity {
public:
    std::string originalFilePath;
    std::string compareFilePath;
    float sampleRate;

    AudioFile<float> originalFile;
    AudioFile<float> compareFile;

    Similarity(std::string originalFilePath, std::string compareFilePath, float sampleRate = 44100.f);
    bool loadAudioFiles();

    float zcrSimilarity();
    float rhythmSimilarity();
    float chromaSimilarity();
    float energyEnvelopeSimilarity();
    float spectralContrastSimilarity();
    float perceptualSimilarity();

    float computeZCR(std::vector<float> buffer);
};

#endif
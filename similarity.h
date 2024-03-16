#ifndef _SIMILARITY_H_
#define _SIMILARITY_H_

#include <string>

#include "AudioFile.h"

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

    float pearsonCorrelation(std::vector<float> x, std::vector<float> y);
    std::vector<float> energyDifference(std::vector<float> buffer);

    std::vector<std::vector<float>> computeSpectrogram(const std::vector<float>& signal, int windowSize, int hopSize);
    std::vector<float> computeSpectralContrast(const std::vector<std::vector<float>>& spectrogram, int numBands);
};

#endif
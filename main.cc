#include <iostream>
#include <string>

#include "similarity.h"

using namespace std;

int main(int argc, char **argv) {
    std::string oboeFilePath = "resources/Oboe.wav";
    std::string trumpetFilePath = "resources/Trumpet.wav";
    float sampleRate = 44100.0;

    Similarity similarity(oboeFilePath, trumpetFilePath, sampleRate);

    float zcrSimilarity = similarity.zcrSimilarity();
    std::cout << "Zero Crossing Rate Similarity: " << zcrSimilarity << std::endl;

    float rhythmSimilarity = similarity.rhythmSimilarity();
    std::cout << "Rhythm Similarity: " << rhythmSimilarity << std::endl;

    float spectralContrastSimilarity = similarity.spectralContrastSimilarity();
    std::cout << "Spectral Contrast Similarity: " << spectralContrastSimilarity << std::endl;

    return 0;
}
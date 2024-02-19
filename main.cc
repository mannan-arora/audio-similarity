#include "similarity.h"
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char **argv) {
    std::string oboeFilePath = "resources/Oboe.wav";
    std::string trumpetFilePath = "resources/Trumpet.wav";
    float sampleRate = 44100.0;

    Similarity similarity(oboeFilePath, trumpetFilePath, sampleRate);

    float zcrSimilarity = similarity.zcrSimilarity();

    std::cout << "Zero Crossing Rate Similarity: " << zcrSimilarity << std::endl;

    return 0;
}
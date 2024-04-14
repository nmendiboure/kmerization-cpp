//
// Created by nmendiboure on 19/12/23.
//

#include <fstream>
#include <sstream>
#include <algorithm>
#include <utility>


#include "genome.h"

// Default constructor implementation
Genome::Genome() {
    this->name_ = "";
    this->description_ = "";
    this->binSize_ = 0;
    this->chrSizes_ = std::map<int, int>();
    this->chrSequences_ = std::map<int, std::pair<std::string, std::string>>();
    this->chrBins_ = Matrix<int>();
    this->idxToBins_ = std::unordered_map<int, std::string>();
    this->binsToIdx_ = std::unordered_map<std::string, int>();
}

Genome::Genome(std::string name) {
    this->name_ = std::move(name);
    this->description_ = "";
    this->binSize_ = 0;
    this->chrSizes_ = std::map<int, int>();
    this->chrSequences_ = std::map<int, std::pair<std::string, std::string>>();
    this->chrBins_ = Matrix<int>();
    this->idxToBins_ = std::unordered_map<int, std::string>();
    this->binsToIdx_ = std::unordered_map<std::string, int>();
}

// Constructor implementation
Genome::Genome(std::string name, std::string  description){
    this->name_ = std::move(name);
    this->description_ = std::move(description);
    this->binSize_ = 0;
    this->chrSizes_ = std::map<int, int>();
    this->chrSequences_ = std::map<int, std::pair<std::string, std::string>>();
    this->chrBins_ = Matrix<int>();
    this->idxToBins_ = std::unordered_map<int, std::string>();
    this->binsToIdx_ = std::unordered_map<std::string, int>();
}

Genome::Genome(std::string name, std::string  description, int binSize){
    this->name_ = std::move(name);
    this->description_ = std::move(description);
    this->binSize_ = binSize;
    this->chrSizes_ = std::map<int, int>();
    this->chrSequences_ = std::map<int, std::pair<std::string, std::string>>();
    this->chrBins_ = Matrix<int>();
    this->idxToBins_ = std::unordered_map<int, std::string>();
    this->binsToIdx_ = std::unordered_map<std::string, int>();
}


int Genome::getGenomeSize() const {
    int genomeSize = 0;
    for (const auto& chr : this->chrSizes_) {
        genomeSize += chr.second;
    }
    return genomeSize;
}




std::string Genome::getChrSequence(int chrId, int strand) const {
    if (strand == 0) return this->chrSequences_.at(chrId).first;
    else return this->chrSequences_.at(chrId).second;
}

void Genome::readFile(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");
    }

    std::string line, chr_id, sequence;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!chr_id.empty()) {
                int id = std::stoi(chr_id.substr(3));
                processSequence(id, sequence);
            }
            std::size_t spacePos = line.find(' ');
            chr_id = line.substr(1, spacePos - 1);
        } else {
            sequence += line;
        }
    }

    if (!chr_id.empty()) {
        int id = std::stoi(chr_id.substr(3));
        processSequence(id, sequence);
    }
}

void Genome::setBinSize(int binSize) {
    this->binSize_ = binSize;
}

void Genome::splitChrInBins() {
    int binCounter = 0;
    for (const auto& chr : this->chrSizes_) {
        int chrId = chr.first;
        int size = chr.second;
        int nbBins = size / binSize_ + 1;
        for (int bin = 0; bin < nbBins; ++bin) {
            int start = bin * binSize_;
            this->chrBins_.push_back(std::vector<int>{chrId, start});
            std::string binName = "chr" + std::to_string(chrId) + ":" + std::to_string(start);
            this->idxToBins_[binCounter] = binName;
            this->binsToIdx_[binName] = binCounter;
            binCounter++;
        }
    }
}


std::string Genome::reverseComplement(const std::string& seq) {
    std::string reversed(seq.rbegin(), seq.rend());
    std::transform(reversed.begin(), reversed.end(), reversed.begin(),
                   [](char c) {
                       switch (c) {
                           case 'a': return 't';
                           case 't': return 'a';
                           case 'g': return 'c';
                           case 'c': return 'g';
                           default: return c;
                       }
                   });
    return reversed;
}


void Genome::processSequence(int id, std::string& seq) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::tolower);
    this->chrSequences_[id] = std::make_pair(seq, reverseComplement(seq));
    this->chrSizes_[id] = seq.length();
    seq.clear();
}

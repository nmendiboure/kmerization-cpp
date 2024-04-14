//
// Created by nicolas on 21/12/23.
//


#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <map>
#include <cmath>
#include <omp.h>
#include "genome.h"
#include "core/component/matrix.h"
#include "core/component/index.h"

#define THREADS 16


// header
std::vector<std::string> splitSequence(const std::string& sequence, int kmerSize);
std::vector<std::string> generateKmers(int KmerSize, std::string alphabet = "atcg");
std::vector<int> IndexKmer(const std::string& sequence, const std::string& kmer);

std::vector<double> getHicFreq(const std::string& filename);
Matrix<double> extractSubMatrix(const Matrix<double>& fullKmersMat, const std::vector<std::string>& fulKmersList, const std::vector<std::string>& subKmersList);
std::tuple<Index, Matrix<double>> computeAllKmers(const Genome& genome, int binSize, const std::vector<std::string>& kmers, const Matrix<int>& chrBins, const std::vector<double>& hicVec);


std::vector<std::string> splitSequence(const std::string& sequence, int kmerSize) {
    int numKmers = static_cast<int>(sequence.size()) - kmerSize + 1;
    std::vector<std::string> kmerList(numKmers);
    for (int i = 0; i < numKmers; ++i) {
        kmerList[i] = sequence.substr(i, kmerSize);
    }
    return kmerList;
}

std::vector<std::string> generateKmers(int KmerSize, std::string alphabet) {
    std::vector<std::string> kmers;
    int nbKmers = pow(static_cast<int>(alphabet.size()), KmerSize);
    kmers.reserve(nbKmers);
    for (int i = 0; i < nbKmers; ++i) {
        std::string kmer;
        int j = i;
        while (j > 0) {
            kmer += alphabet[j % alphabet.size()];
            j /= alphabet.size();
        }
        while (kmer.size() < KmerSize) {
            kmer += alphabet[0];
        }
        kmers.emplace_back(kmer);
    }

    return kmers;
}

std::vector<int> IndexKmer(const std::string& sequence, const std::string& kmer) {
    std::vector<int> positions;
    std::string::size_type pos = 0;

    while (true) {
        pos = sequence.find(kmer, pos);
        if (pos == std::string::npos) {
            // No more occurrences found, exit the loop
            break;
        }
        positions.emplace_back(static_cast<int>(pos));
        pos++; // Move to the next character for overlapping matches
    }

    return positions;
}


std::tuple<Index, Matrix<double>> computeAllKmers(
        const Genome& genome,
        int binSize,
        const std::vector<std::string>& kmers,
        const Matrix<int>& chrBins,
        const std::vector<double>& hicVec
        ){

    // Precompute all necessary sequences
    std::vector<std::pair<std::string, std::string>> chrSequences;
    std::vector<std::string> sortedKmers = kmers;
    std::sort(sortedKmers.begin(), sortedKmers.end());


    size_t nbChr = genome.getChrCount();
    for (int chr = 0; chr < nbChr; ++chr) {
        int chrId = chr + 1;
        chrSequences.emplace_back(genome.getChrSequence(chrId, 0),
                                  genome.getChrSequence(chrId, 1));
    }

    // Create a vector of unordered_maps for each thread
    std::vector<std::unordered_map<std::string, Matrix<int>>> local_maps(THREADS);

    Matrix<double> kmersMat1(static_cast<int>(kmers.size()), static_cast<int>(chrBins.getHeight()),0.0);

    int bin;
    #pragma omp parallel for schedule(dynamic) num_threads(THREADS)
    for (int k = 0; k < sortedKmers.size(); ++k) {
        const auto& kmer = sortedKmers[k];
        int threadId = omp_get_thread_num();
        auto& local_map = local_maps[threadId];
        Matrix<int> indexMat;

        for (int chr = 0; chr < nbChr; ++chr) {
            int chrId = chr + 1;
            const auto& watsonSeq = chrSequences[chr].first;
            const auto& crickSeq = chrSequences[chr].second;

            auto watsonPositions = IndexKmer(watsonSeq, kmer);
            auto crickPositions = IndexKmer(crickSeq, kmer);

            indexMat.reserve(watsonPositions.size() + crickPositions.size(), 4);
            for (const auto& w_pos : watsonPositions) {
                bin = w_pos / binSize * binSize;
                indexMat.push_back({chrId, 1, w_pos, bin});
            }
            for (const auto& c_pos : crickPositions) {
                int reversed_c_pos = genome.getChrSize(chrId) - c_pos - kmer.size();
                bin =  reversed_c_pos / binSize * binSize;
                indexMat.push_back({chrId, -1, reversed_c_pos, bin});
            }
        }

        for (int i = 0; i < chrBins.getHeight(); ++i) {
            const auto& chrBin = chrBins(i);
            int count = 0;
            for (int j = 0; j < indexMat.getHeight(); ++j) {
                if (indexMat(j, 0) == chrBin[0] && indexMat(j, 3) == chrBin[1]) {
                    count++;
                }
            }

            kmersMat1(k, i) = count;
        }

        local_map[kmer] = std::move(indexMat);
    }

    Index KmersIndex;
    for (const auto& local_map : local_maps) {
        for (const auto& kv : local_map) {
            KmersIndex.insert(kv.first, kv.second.duplicate());
        }
    }

    Matrix<double> kmersMat2 = kmersMat1.multiply(hicVec);
    Matrix<double> kmersMat3 = kmersMat2.divide(kmersMat2.sum());
    Matrix<double> kmersMat = kmersMat3.cumuSum(1);

    return {std::move(KmersIndex), std::move(kmersMat)};
}


std::vector<double> getHicFreq(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("File not found or could not be opened.");
    }

    std::string line;
    std::vector<double> frequencies;
    frequencies.reserve(10000);  // Reserve some capacity if you expect a large file

    // Read headers and find the index of the 'frequencies' column
    if (!std::getline(file, line)) {
        throw std::runtime_error("File is empty or header is missing.");
    }
    std::istringstream headerStream(line);
    std::string cell;
    int frequencyIndex = -1;
    int columnIndex = 0;
    while (std::getline(headerStream, cell, '\t')) {
        if (cell == "frequencies") {
            frequencyIndex = columnIndex;
            break;
        }
        columnIndex++;
    }
    if (frequencyIndex == -1) {
        throw std::runtime_error("Column 'frequencies' not found in the file.");
    }

    // Read the frequencies column
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        for (int i = 0; i <= frequencyIndex; ++i) {
            if (!std::getline(lineStream, cell, '\t')) {
                throw std::runtime_error("Error reading frequency value from line.");
            }
        }
        frequencies.emplace_back(std::stod(cell));
    }

    return frequencies;
}


Matrix<double> extractSubMatrix(
        const Matrix<double>& fullKmersMat,
        const std::vector<std::string>& fulKmersList,
        const std::vector<std::string>& subKmersList)
        {

    Matrix<double> subMat;
    std::vector<std::string> sortedFulKmersList = fulKmersList;
    std::sort(sortedFulKmersList.begin(), sortedFulKmersList.end());

    // Iterate over subKmersList and find each kmer in sortedFulKmersList
    for (const auto& subKmer : subKmersList) {
        auto it = std::lower_bound(sortedFulKmersList.begin(), sortedFulKmersList.end(), subKmer);
        if (it != sortedFulKmersList.end() && *it == subKmer) {
            // If found, get the index and add the corresponding row to subMat
            int index = std::distance(sortedFulKmersList.begin(), it);
            subMat.push_back(fullKmersMat(index));
        }
    }

    return subMat;
}

int main(int argc, char* argv[])
{

    // Check for the --clear argument
    bool clearFlag = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--clear") {
            clearFlag = true;
            break;
        }
    }

    // paths
    std::string genomePath = "../data/inputs/S288c.fa";
    std::string hicPath = "../data/inputs/AD403_20kb.tsv";
    std::string kmersIndexPath = "../data/inputs/kmersIndex20kb.bin";
    std::string fullMatPath = "../data/inputs/fullMat20kb.bin";
    std::string LyMatPath = "../data/inputs/LY_Mat20kb.bin";
    std::string LyIndexPath = "../data/inputs/LY_Index20kb.bin";

    // Delete existing files if --clear flag is set
    if (clearFlag) {
        std::cout << "Clearing existing binary files..." << std::endl;
        std::filesystem::remove(kmersIndexPath);
        std::filesystem::remove(fullMatPath);
        std::filesystem::remove(LyMatPath);
        std::filesystem::remove(LyIndexPath);
    }


    int kmerSize = 8;
    int binSize = 20000;
    std::string LySeq = "atgactaacgaaaaggtctggatagagaagttggataatccaactctttcagtgttaccacatgactttttacgcccacaacaagaaccttatacgaaacaagctacatattcgttacagctacctcagctcgatgtgcctcatgatagtttttctaacaaatacgctgtcgctttgagtgtatgggctgcattgatatatagagtaaccggtgacgatgatattgttctttatattgcgaataacaaaatcttaagattcaatattcaaccaacgtggtcatttaatgagctgtattctacaattaacaatgagttgaacaagctcaattctattgaggccaatttttcctttgacgagctagctgaaaaaattcaaagttgccaagatctggaaaggacccctcagttgttccgtttggcctttttggaaaaccaagatttcaaattagacgagttcaagcatcatttagtggactttgctttgaatttggataccagtaataatgcgcatgttttgaacttaatttataacagcttactgtattcgaatgaaagagtaaccattgttgcggaccaatttactcaatatttgactgctgcgctaagcgatccatccaattgcataactaaaatctctctgatcaccgcatcatccaaggatagtttacctgatccaactaagaacttgggctggtgcgatttcgtggggtgtattcacgacattttccaggacaatgctgaagccttcccagagagaacctgtgttgtggagactccaacactaaattccgacaagtcccgttctttcacttatcgcgacatcaaccgcacttctaacatagttgcccattatttgattaaaacaggtatcaaaagaggtgatgtagtgatgatctattcttctaggggtgtggatttgatggtatgtgtgatgggtgtcttgaaagccggcgcaaccttttcagttatcgaccctgcatatcccccagccagacaaaccatttacttaggtgttgctaaaccacgtgggttgattgttattagagctgctggacaattggatcaactagtagaagattacatcaatgatgaattggagattgtttcaagaatcaattccatcgctattcaagaaaatggtaccattgaaggtggcaaattggacaatggcgaggatgttttggctccatatgatcactacaaagacaccagaacaggtgttgtagttggaccagattccaacccaaccctatctttcacatctggttccgaaggtattcctaagggtgttcttggtagacatttttccttggcttattatttcaattggatgtccaaaaggttcaacttaacagaaaatgataaattcacaatgctgagcggtattgcacatgatccaattcaaagagatatgtttacaccattatttttaggtgcccaattgtatgtccctactcaagatgatattggtacaccgggccgtttagcggaatggatgagtaagtatggttgcacagttacccatttaacacctgccatgggtcaattacttactgcccaagctactacaccattccctaagttacatcatgcgttctttgtgggtgacattttaacaaaacgtgattgtctgaggttacaaaccttggcagaaaattgccgtattgttaatatgtacggtaccactgaaacacagcgtgcagtttcttatttcgaagttaaatcaaaaaatgacgatccaaactttttgaaaaaattgaaagatgtcatgcctgctggtaaaggtatgttgaacgttcagctactagttgttaacaggaacgatcgtactcaaatatgtggtattggcgaaataggtgagatttatgttcgtgcaggtggtttggccgaaggttatagaggattaccagaattgaataaagaaaaatttgtgaacaactggtttgttgaaaaagatcactggaattatttggataaggataatggtgaaccttggagacaattctggttaggtccaagagatagattgtacagaacgggtgatttaggtcgttatctaccaaacgg";
    Genome myGenome("S288C", "Saccharomyces cerevisiae S288C", binSize);
    myGenome.readFile(genomePath);
    myGenome.splitChrInBins();

    const Matrix<int>& chrBins = myGenome.getChrBinsMat();
    std::vector<std::string> LySeqList = splitSequence(LySeq, kmerSize);
    std::vector<std::string> kmersList = generateKmers(kmerSize, "atcg");

//    kmersList.erase(kmersList.begin() + 10000, kmersList.end());

    std::cout << "Import Hi-C data" << std::endl;
    std::vector<double> hicData = getHicFreq(hicPath);

    Index kmersIndex, LyIndex;
    Matrix<double> fullMatrix, LyMat;

    if (!std::filesystem::exists(kmersIndexPath)) {
        std::cout << "Indexing all kmers" << std::endl;
        std::tie(kmersIndex, fullMatrix) = computeAllKmers(myGenome, binSize, kmersList, chrBins, hicData);
        std::cout << "Dumping kmers index to binary file" << std::endl;
        kmersIndex.dumpToProto(kmersIndexPath);
        std::cout << "Dumping full core to binary file" << std::endl;
        fullMatrix.dumpToProto(fullMatPath);
    } else {
        std::cout << "Loading kmers index from binary file" << std::endl;
        kmersIndex = Index::loadFromProto(kmersIndexPath);
        std::cout << "Loading full core from binary file" << std::endl;
        fullMatrix = Matrix<double>::loadFromProto(fullMatPath);
    }

    if (!std::filesystem::exists(LyMatPath)) {
        std::cout << "Extracting LY kmers" << std::endl;
        LyIndex = kmersIndex.subSet(LySeqList);
        LyMat = extractSubMatrix(fullMatrix, kmersList, LySeqList);
        std::cout << "Dumping LY core to binary file" << std::endl;
        LyMat.dumpToProto(LyMatPath);
        std::cout << "Dumping LY index to binary file" << std::endl;
        LyIndex.dumpToProto(LyIndexPath);
    } else {
        std::cout << "Loading LY core from binary file" << std::endl;
        LyMat = Matrix<double>::loadFromProto(LyMatPath);
        std::cout << "Loading LY index from binary file" << std::endl;
        LyIndex= Index::loadFromProto(LyIndexPath);
    }

    return 0;
}


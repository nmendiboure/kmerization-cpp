//
// Created by nmendiboure on 19/12/23.
//

#ifndef GENOME_H
#define GENOME_H

#include <string>
#include <unordered_map>
#include <map>
#include "core/component/matrix.h"

class Genome {
public:
    Genome();
    explicit Genome(std::string  name);
    Genome(std::string  name, std::string  description);
    Genome(std::string  name, std::string  description, int binSize);

    Genome(Genome& genome) = delete;

    // Destructor
    virtual ~Genome() = default;

    // Method to read a FASTA file
    void readFile(const std::string& path);

    // Public attributes
    inline int getChrCount() const { return static_cast<int>(chrSizes_.size()); }
    inline std::string getName() const { return name_; }
    inline std::string getDescription() const { return description_; }
    inline int getBinSize() const { return binSize_; }
    inline const Matrix<int>& getChrBinsMat() { return chrBins_; }
    inline int getChrSize(int chrId) const { return chrSizes_.at(chrId); }

    void setBinSize(int binSize);

    int getGenomeSize() const;
    std::string getChrSequence(int chrId, int strand) const;

    // Methods
    void splitChrInBins();
    static std::string reverseComplement(const std::string& seq);

private:
    std::string name_;
    std::string description_;
    int binSize_;

    std::map<int, int> chrSizes_;
    std::map<int, std::pair<std::string, std::string>> chrSequences_;

    Matrix<int> chrBins_;
    std::unordered_map<int, std::string> idxToBins_;
    std::unordered_map<std::string, int> binsToIdx_;

    void processSequence(int id, std::string &seq);


};

#endif // GENOME_H


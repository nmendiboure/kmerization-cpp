//
// Created by nicolas on 23/12/23.
//

#include "index.h"
#include <fstream>
#include <memory>

// default constructor
Index::Index() : data() {}

// constructor from a map of string - core pairs
Index::Index(const std::unordered_map<std::string, std::unique_ptr<Matrix<int>>>& data) {
    for (const auto& pair : data) {
        this->data[pair.first] = std::make_unique<Matrix<int>>(pair.second->duplicate());
    }
}
// move constructor
Index::Index(Index&& other) noexcept : data(std::move(other.data)) {}


// insertion of a new key - value pair
void Index::insert(const std::string& key, Matrix<int>&& value) {
    data[key] = std::make_unique<Matrix<int>>(std::move(value));;
}


// get a subset of the index
Index Index::subSet(const std::vector<std::string>& keys) const {
    Index subIndex;
    for (const auto& key : keys) {
        if (this->contains(key)) {
            Matrix<int> matrixCopy = this->get(key)->duplicate();
            subIndex.insert(key, std::move(matrixCopy));
        }
    }
    return subIndex;
}

// get all the keys
std::vector<std::string> Index::listKeys() const {
    std::vector<std::string> keysList;
    keysList.reserve(data.size());
    for (const auto& pair : data) {
        keysList.push_back(pair.first);
    }

    return keysList;
}

// get all the key - value pairs
std::vector<std::pair<std::string, const Matrix<int>*>> Index::items() const {
    std::vector<std::pair<std::string, const Matrix<int>*>> itemsList;
    itemsList.reserve(data.size());
    for (const auto& pair : data) {
        itemsList.emplace_back(pair.first, pair.second.get());
    }
    return itemsList;
}

// Move assignment operator definition
Index& Index::operator=(Index&& other) noexcept {
    if (this != &other) {
        data = std::move(other.data);
    }
    return *this;
}


void Index::dumpToProto(const std::string& filePath) const {
    protoIndex protoInd;
    IndexToProto(*this, protoInd);
    std::ofstream outFile(filePath, std::ios::binary);
    protoInd.SerializeToOstream(&outFile);
    outFile.close();
}

Index Index::loadFromProto(const std::string& filePath) {
    protoIndex protoInd;
    std::ifstream inFile(filePath, std::ios::binary);
    protoInd.ParseFromIstream(&inFile);
    inFile.close();

    return ProtoToIndex(protoInd);
}

void IndexToProto(const Index& index, protoIndex& protoInd){
    for (const auto& kv : index.items()) {
        protoIndexEntry* entry = protoInd.add_entries();
        entry->set_key(kv.first);
        MatrixToProto(*kv.second, *entry->mutable_value());
    }
}


Index ProtoToIndex(const protoIndex& protoInd){
    Index index;
    for (const auto& entry : protoInd.entries()) {
        index.insert(entry.key(), ProtoToMatrix<int>(entry.value()));
    }

    return index;
}
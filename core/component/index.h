#ifndef INDEX_H
#define INDEX_H

#include <string>
#include <unordered_map>
#include <memory>
#include "matrix.h"
#include "proto/index.pb.h"

class Index {
public:
    Index();
    explicit Index(const std::unordered_map<std::string, std::unique_ptr<Matrix<int>>>& data);
    Index(Index&& other) noexcept;
    Index(const Index& other) = delete;

    virtual ~Index() = default;

    inline Matrix<int>& operator()(const std::string& key) { return *data.at(key); }
    inline const Matrix<int>& operator()(const std::string& key) const { return *data.at(key); }
    inline Matrix<int>* get(const std::string& key) { return data.at(key).get(); }
    inline const Matrix<int>* get(const std::string& key) const { return data.at(key).get(); }

    void insert(const std::string& key, Matrix<int>&& value);
    inline bool contains(const std::string& key) const { return data.find(key) != data.end(); }
    Index subSet(const std::vector<std::string>& keys) const;

    std::vector<std::string> listKeys() const;
    std::vector<std::pair<std::string, const Matrix<int>*>> items() const;

    inline size_t getSize() const { return data.size(); }
    inline bool isEmpty() const { return data.empty(); }
    inline void clear() { data.clear(); }

    Index& operator=(Index&& other) noexcept;

    void dumpToProto(const std::string& filePath) const;
    static Index loadFromProto(const std::string& filePath);

private:
    std::unordered_map<std::string, std::unique_ptr<Matrix<int>>> data;
};

void IndexToProto(const Index& index, protoIndex& protoInd);
Index ProtoToIndex(const protoIndex& protoInd);

#endif // INDEX_H
// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: matrix.proto

#include "matrix.pb.h"

#include <algorithm>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/wire_format_lite.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>

PROTOBUF_PRAGMA_INIT_SEG

namespace _pb = ::PROTOBUF_NAMESPACE_ID;
namespace _pbi = _pb::internal;

PROTOBUF_CONSTEXPR protoMatrix::protoMatrix(
    ::_pbi::ConstantInitialized): _impl_{
    /*decltype(_impl_.data_)*/{}
  , /*decltype(_impl_.height_)*/0
  , /*decltype(_impl_.width_)*/0
  , /*decltype(_impl_._cached_size_)*/{}} {}
struct protoMatrixDefaultTypeInternal {
  PROTOBUF_CONSTEXPR protoMatrixDefaultTypeInternal()
      : _instance(::_pbi::ConstantInitialized{}) {}
  ~protoMatrixDefaultTypeInternal() {}
  union {
    protoMatrix _instance;
  };
};
PROTOBUF_ATTRIBUTE_NO_DESTROY PROTOBUF_CONSTINIT PROTOBUF_ATTRIBUTE_INIT_PRIORITY1 protoMatrixDefaultTypeInternal _protoMatrix_default_instance_;
static ::_pb::Metadata file_level_metadata_matrix_2eproto[1];
static constexpr ::_pb::EnumDescriptor const** file_level_enum_descriptors_matrix_2eproto = nullptr;
static constexpr ::_pb::ServiceDescriptor const** file_level_service_descriptors_matrix_2eproto = nullptr;

const uint32_t TableStruct_matrix_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::protoMatrix, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  ~0u,  // no _inlined_string_donated_
  PROTOBUF_FIELD_OFFSET(::protoMatrix, _impl_.height_),
  PROTOBUF_FIELD_OFFSET(::protoMatrix, _impl_.width_),
  PROTOBUF_FIELD_OFFSET(::protoMatrix, _impl_.data_),
};
static const ::_pbi::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, -1, -1, sizeof(::protoMatrix)},
};

static const ::_pb::Message* const file_default_instances[] = {
  &::_protoMatrix_default_instance_._instance,
};

const char descriptor_table_protodef_matrix_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\014matrix.proto\":\n\013protoMatrix\022\016\n\006height\030"
  "\001 \001(\005\022\r\n\005width\030\002 \001(\005\022\014\n\004data\030\003 \003(\001b\006prot"
  "o3"
  ;
static ::_pbi::once_flag descriptor_table_matrix_2eproto_once;
const ::_pbi::DescriptorTable descriptor_table_matrix_2eproto = {
    false, false, 82, descriptor_table_protodef_matrix_2eproto,
    "matrix.proto",
    &descriptor_table_matrix_2eproto_once, nullptr, 0, 1,
    schemas, file_default_instances, TableStruct_matrix_2eproto::offsets,
    file_level_metadata_matrix_2eproto, file_level_enum_descriptors_matrix_2eproto,
    file_level_service_descriptors_matrix_2eproto,
};
PROTOBUF_ATTRIBUTE_WEAK const ::_pbi::DescriptorTable* descriptor_table_matrix_2eproto_getter() {
  return &descriptor_table_matrix_2eproto;
}

// Force running AddDescriptors() at dynamic initialization time.
PROTOBUF_ATTRIBUTE_INIT_PRIORITY2 static ::_pbi::AddDescriptorsRunner dynamic_init_dummy_matrix_2eproto(&descriptor_table_matrix_2eproto);

// ===================================================================

class protoMatrix::_Internal {
 public:
};

protoMatrix::protoMatrix(::PROTOBUF_NAMESPACE_ID::Arena* arena,
                         bool is_message_owned)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena, is_message_owned) {
  SharedCtor(arena, is_message_owned);
  // @@protoc_insertion_point(arena_constructor:protoMatrix)
}
protoMatrix::protoMatrix(const protoMatrix& from)
  : ::PROTOBUF_NAMESPACE_ID::Message() {
  protoMatrix* const _this = this; (void)_this;
  new (&_impl_) Impl_{
      decltype(_impl_.data_){from._impl_.data_}
    , decltype(_impl_.height_){}
    , decltype(_impl_.width_){}
    , /*decltype(_impl_._cached_size_)*/{}};

  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&_impl_.height_, &from._impl_.height_,
    static_cast<size_t>(reinterpret_cast<char*>(&_impl_.width_) -
    reinterpret_cast<char*>(&_impl_.height_)) + sizeof(_impl_.width_));
  // @@protoc_insertion_point(copy_constructor:protoMatrix)
}

inline void protoMatrix::SharedCtor(
    ::_pb::Arena* arena, bool is_message_owned) {
  (void)arena;
  (void)is_message_owned;
  new (&_impl_) Impl_{
      decltype(_impl_.data_){arena}
    , decltype(_impl_.height_){0}
    , decltype(_impl_.width_){0}
    , /*decltype(_impl_._cached_size_)*/{}
  };
}

protoMatrix::~protoMatrix() {
  // @@protoc_insertion_point(destructor:protoMatrix)
  if (auto *arena = _internal_metadata_.DeleteReturnArena<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>()) {
  (void)arena;
    return;
  }
  SharedDtor();
}

inline void protoMatrix::SharedDtor() {
  GOOGLE_DCHECK(GetArenaForAllocation() == nullptr);
  _impl_.data_.~RepeatedField();
}

void protoMatrix::SetCachedSize(int size) const {
  _impl_._cached_size_.Set(size);
}

void protoMatrix::Clear() {
// @@protoc_insertion_point(message_clear_start:protoMatrix)
  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  _impl_.data_.Clear();
  ::memset(&_impl_.height_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&_impl_.width_) -
      reinterpret_cast<char*>(&_impl_.height_)) + sizeof(_impl_.width_));
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* protoMatrix::_InternalParse(const char* ptr, ::_pbi::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  while (!ctx->Done(&ptr)) {
    uint32_t tag;
    ptr = ::_pbi::ReadTag(ptr, &tag);
    switch (tag >> 3) {
      // int32 height = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 8)) {
          _impl_.height_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // int32 width = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 16)) {
          _impl_.width_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint32(&ptr);
          CHK_(ptr);
        } else
          goto handle_unusual;
        continue;
      // repeated double data = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<uint8_t>(tag) == 26)) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_data(), ptr, ctx);
          CHK_(ptr);
        } else if (static_cast<uint8_t>(tag) == 25) {
          _internal_add_data(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
          ptr += sizeof(double);
        } else
          goto handle_unusual;
        continue;
      default:
        goto handle_unusual;
    }  // switch
  handle_unusual:
    if ((tag == 0) || ((tag & 7) == 4)) {
      CHK_(ptr);
      ctx->SetLastTag(tag);
      goto message_done;
    }
    ptr = UnknownFieldParse(
        tag,
        _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
        ptr, ctx);
    CHK_(ptr != nullptr);
  }  // while
message_done:
  return ptr;
failure:
  ptr = nullptr;
  goto message_done;
#undef CHK_
}

uint8_t* protoMatrix::_InternalSerialize(
    uint8_t* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:protoMatrix)
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  // int32 height = 1;
  if (this->_internal_height() != 0) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteInt32ToArray(1, this->_internal_height(), target);
  }

  // int32 width = 2;
  if (this->_internal_width() != 0) {
    target = stream->EnsureSpace(target);
    target = ::_pbi::WireFormatLite::WriteInt32ToArray(2, this->_internal_width(), target);
  }

  // repeated double data = 3;
  if (this->_internal_data_size() > 0) {
    target = stream->WriteFixedPacked(3, _internal_data(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::_pbi::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:protoMatrix)
  return target;
}

size_t protoMatrix::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:protoMatrix)
  size_t total_size = 0;

  uint32_t cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated double data = 3;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_data_size());
    size_t data_size = 8UL * count;
    if (data_size > 0) {
      total_size += 1 +
        ::_pbi::WireFormatLite::Int32Size(static_cast<int32_t>(data_size));
    }
    total_size += data_size;
  }

  // int32 height = 1;
  if (this->_internal_height() != 0) {
    total_size += ::_pbi::WireFormatLite::Int32SizePlusOne(this->_internal_height());
  }

  // int32 width = 2;
  if (this->_internal_width() != 0) {
    total_size += ::_pbi::WireFormatLite::Int32SizePlusOne(this->_internal_width());
  }

  return MaybeComputeUnknownFieldsSize(total_size, &_impl_._cached_size_);
}

const ::PROTOBUF_NAMESPACE_ID::Message::ClassData protoMatrix::_class_data_ = {
    ::PROTOBUF_NAMESPACE_ID::Message::CopyWithSourceCheck,
    protoMatrix::MergeImpl
};
const ::PROTOBUF_NAMESPACE_ID::Message::ClassData*protoMatrix::GetClassData() const { return &_class_data_; }


void protoMatrix::MergeImpl(::PROTOBUF_NAMESPACE_ID::Message& to_msg, const ::PROTOBUF_NAMESPACE_ID::Message& from_msg) {
  auto* const _this = static_cast<protoMatrix*>(&to_msg);
  auto& from = static_cast<const protoMatrix&>(from_msg);
  // @@protoc_insertion_point(class_specific_merge_from_start:protoMatrix)
  GOOGLE_DCHECK_NE(&from, _this);
  uint32_t cached_has_bits = 0;
  (void) cached_has_bits;

  _this->_impl_.data_.MergeFrom(from._impl_.data_);
  if (from._internal_height() != 0) {
    _this->_internal_set_height(from._internal_height());
  }
  if (from._internal_width() != 0) {
    _this->_internal_set_width(from._internal_width());
  }
  _this->_internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
}

void protoMatrix::CopyFrom(const protoMatrix& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:protoMatrix)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool protoMatrix::IsInitialized() const {
  return true;
}

void protoMatrix::InternalSwap(protoMatrix* other) {
  using std::swap;
  _internal_metadata_.InternalSwap(&other->_internal_metadata_);
  _impl_.data_.InternalSwap(&other->_impl_.data_);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(protoMatrix, _impl_.width_)
      + sizeof(protoMatrix::_impl_.width_)
      - PROTOBUF_FIELD_OFFSET(protoMatrix, _impl_.height_)>(
          reinterpret_cast<char*>(&_impl_.height_),
          reinterpret_cast<char*>(&other->_impl_.height_));
}

::PROTOBUF_NAMESPACE_ID::Metadata protoMatrix::GetMetadata() const {
  return ::_pbi::AssignDescriptors(
      &descriptor_table_matrix_2eproto_getter, &descriptor_table_matrix_2eproto_once,
      file_level_metadata_matrix_2eproto[0]);
}

// @@protoc_insertion_point(namespace_scope)
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::protoMatrix*
Arena::CreateMaybeMessage< ::protoMatrix >(Arena* arena) {
  return Arena::CreateMessageInternal< ::protoMatrix >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>

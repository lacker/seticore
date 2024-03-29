// Generated by Cap'n Proto compiler, DO NOT EDIT
// source: stamp.capnp

#pragma once

#include <capnp/generated-header-support.h>
#include <kj/windows-sanity.h>

#if CAPNP_VERSION != 11000
#error "Version mismatch between generated code and library headers.  You must use the same version of the Cap'n Proto compiler and library."
#endif

#include "hit.capnp.h"

CAPNP_BEGIN_HEADER

namespace capnp {
namespace schemas {

CAPNP_DECLARE_SCHEMA(ce98ac82349b1161);

}  // namespace schemas
}  // namespace capnp


struct Stamp {
  Stamp() = delete;

  class Reader;
  class Builder;
  class Pipeline;

  struct _capnpPrivate {
    CAPNP_DECLARE_STRUCT_HEADER(ce98ac82349b1161, 11, 5)
    #if !CAPNP_LITE
    static constexpr ::capnp::_::RawBrandedSchema const* brand() { return &schema->defaultBrand; }
    #endif  // !CAPNP_LITE
  };
};

// =======================================================================================

class Stamp::Reader {
public:
  typedef Stamp Reads;

  Reader() = default;
  inline explicit Reader(::capnp::_::StructReader base): _reader(base) {}

  inline ::capnp::MessageSize totalSize() const {
    return _reader.totalSize().asPublic();
  }

#if !CAPNP_LITE
  inline ::kj::StringTree toString() const {
    return ::capnp::_::structString(_reader, *_capnpPrivate::brand());
  }
#endif  // !CAPNP_LITE

  inline bool hasSourceName() const;
  inline  ::capnp::Text::Reader getSourceName() const;

  inline double getRa() const;

  inline double getDec() const;

  inline double getFch1() const;

  inline double getFoff() const;

  inline double getTstart() const;

  inline double getTsamp() const;

  inline  ::int32_t getTelescopeId() const;

  inline  ::int32_t getNumTimesteps() const;

  inline  ::int32_t getNumChannels() const;

  inline  ::int32_t getNumPolarizations() const;

  inline  ::int32_t getNumAntennas() const;

  inline bool hasData() const;
  inline  ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>::Reader getData() const;

  inline bool hasSeticoreVersion() const;
  inline  ::capnp::Text::Reader getSeticoreVersion() const;

  inline  ::int32_t getCoarseChannel() const;

  inline  ::int32_t getFftSize() const;

  inline  ::int32_t getStartChannel() const;

  inline bool hasSignal() const;
  inline  ::Signal::Reader getSignal() const;

  inline  ::int32_t getSchan() const;

  inline bool hasObsid() const;
  inline  ::capnp::Text::Reader getObsid() const;

private:
  ::capnp::_::StructReader _reader;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::ToDynamic_;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::_::PointerHelpers;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::List;
  friend class ::capnp::MessageBuilder;
  friend class ::capnp::Orphanage;
};

class Stamp::Builder {
public:
  typedef Stamp Builds;

  Builder() = delete;  // Deleted to discourage incorrect usage.
                       // You can explicitly initialize to nullptr instead.
  inline Builder(decltype(nullptr)) {}
  inline explicit Builder(::capnp::_::StructBuilder base): _builder(base) {}
  inline operator Reader() const { return Reader(_builder.asReader()); }
  inline Reader asReader() const { return *this; }

  inline ::capnp::MessageSize totalSize() const { return asReader().totalSize(); }
#if !CAPNP_LITE
  inline ::kj::StringTree toString() const { return asReader().toString(); }
#endif  // !CAPNP_LITE

  inline bool hasSourceName();
  inline  ::capnp::Text::Builder getSourceName();
  inline void setSourceName( ::capnp::Text::Reader value);
  inline  ::capnp::Text::Builder initSourceName(unsigned int size);
  inline void adoptSourceName(::capnp::Orphan< ::capnp::Text>&& value);
  inline ::capnp::Orphan< ::capnp::Text> disownSourceName();

  inline double getRa();
  inline void setRa(double value);

  inline double getDec();
  inline void setDec(double value);

  inline double getFch1();
  inline void setFch1(double value);

  inline double getFoff();
  inline void setFoff(double value);

  inline double getTstart();
  inline void setTstart(double value);

  inline double getTsamp();
  inline void setTsamp(double value);

  inline  ::int32_t getTelescopeId();
  inline void setTelescopeId( ::int32_t value);

  inline  ::int32_t getNumTimesteps();
  inline void setNumTimesteps( ::int32_t value);

  inline  ::int32_t getNumChannels();
  inline void setNumChannels( ::int32_t value);

  inline  ::int32_t getNumPolarizations();
  inline void setNumPolarizations( ::int32_t value);

  inline  ::int32_t getNumAntennas();
  inline void setNumAntennas( ::int32_t value);

  inline bool hasData();
  inline  ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>::Builder getData();
  inline void setData( ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>::Reader value);
  inline void setData(::kj::ArrayPtr<const float> value);
  inline  ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>::Builder initData(unsigned int size);
  inline void adoptData(::capnp::Orphan< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>>&& value);
  inline ::capnp::Orphan< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>> disownData();

  inline bool hasSeticoreVersion();
  inline  ::capnp::Text::Builder getSeticoreVersion();
  inline void setSeticoreVersion( ::capnp::Text::Reader value);
  inline  ::capnp::Text::Builder initSeticoreVersion(unsigned int size);
  inline void adoptSeticoreVersion(::capnp::Orphan< ::capnp::Text>&& value);
  inline ::capnp::Orphan< ::capnp::Text> disownSeticoreVersion();

  inline  ::int32_t getCoarseChannel();
  inline void setCoarseChannel( ::int32_t value);

  inline  ::int32_t getFftSize();
  inline void setFftSize( ::int32_t value);

  inline  ::int32_t getStartChannel();
  inline void setStartChannel( ::int32_t value);

  inline bool hasSignal();
  inline  ::Signal::Builder getSignal();
  inline void setSignal( ::Signal::Reader value);
  inline  ::Signal::Builder initSignal();
  inline void adoptSignal(::capnp::Orphan< ::Signal>&& value);
  inline ::capnp::Orphan< ::Signal> disownSignal();

  inline  ::int32_t getSchan();
  inline void setSchan( ::int32_t value);

  inline bool hasObsid();
  inline  ::capnp::Text::Builder getObsid();
  inline void setObsid( ::capnp::Text::Reader value);
  inline  ::capnp::Text::Builder initObsid(unsigned int size);
  inline void adoptObsid(::capnp::Orphan< ::capnp::Text>&& value);
  inline ::capnp::Orphan< ::capnp::Text> disownObsid();

private:
  ::capnp::_::StructBuilder _builder;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::ToDynamic_;
  friend class ::capnp::Orphanage;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::_::PointerHelpers;
};

#if !CAPNP_LITE
class Stamp::Pipeline {
public:
  typedef Stamp Pipelines;

  inline Pipeline(decltype(nullptr)): _typeless(nullptr) {}
  inline explicit Pipeline(::capnp::AnyPointer::Pipeline&& typeless)
      : _typeless(kj::mv(typeless)) {}

  inline  ::Signal::Pipeline getSignal();
private:
  ::capnp::AnyPointer::Pipeline _typeless;
  friend class ::capnp::PipelineHook;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::ToDynamic_;
};
#endif  // !CAPNP_LITE

// =======================================================================================

inline bool Stamp::Reader::hasSourceName() const {
  return !_reader.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS).isNull();
}
inline bool Stamp::Builder::hasSourceName() {
  return !_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS).isNull();
}
inline  ::capnp::Text::Reader Stamp::Reader::getSourceName() const {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::get(_reader.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS));
}
inline  ::capnp::Text::Builder Stamp::Builder::getSourceName() {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::get(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS));
}
inline void Stamp::Builder::setSourceName( ::capnp::Text::Reader value) {
  ::capnp::_::PointerHelpers< ::capnp::Text>::set(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS), value);
}
inline  ::capnp::Text::Builder Stamp::Builder::initSourceName(unsigned int size) {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::init(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS), size);
}
inline void Stamp::Builder::adoptSourceName(
    ::capnp::Orphan< ::capnp::Text>&& value) {
  ::capnp::_::PointerHelpers< ::capnp::Text>::adopt(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS), kj::mv(value));
}
inline ::capnp::Orphan< ::capnp::Text> Stamp::Builder::disownSourceName() {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::disown(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS));
}

inline double Stamp::Reader::getRa() const {
  return _reader.getDataField<double>(
      ::capnp::bounded<0>() * ::capnp::ELEMENTS);
}

inline double Stamp::Builder::getRa() {
  return _builder.getDataField<double>(
      ::capnp::bounded<0>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setRa(double value) {
  _builder.setDataField<double>(
      ::capnp::bounded<0>() * ::capnp::ELEMENTS, value);
}

inline double Stamp::Reader::getDec() const {
  return _reader.getDataField<double>(
      ::capnp::bounded<1>() * ::capnp::ELEMENTS);
}

inline double Stamp::Builder::getDec() {
  return _builder.getDataField<double>(
      ::capnp::bounded<1>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setDec(double value) {
  _builder.setDataField<double>(
      ::capnp::bounded<1>() * ::capnp::ELEMENTS, value);
}

inline double Stamp::Reader::getFch1() const {
  return _reader.getDataField<double>(
      ::capnp::bounded<2>() * ::capnp::ELEMENTS);
}

inline double Stamp::Builder::getFch1() {
  return _builder.getDataField<double>(
      ::capnp::bounded<2>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setFch1(double value) {
  _builder.setDataField<double>(
      ::capnp::bounded<2>() * ::capnp::ELEMENTS, value);
}

inline double Stamp::Reader::getFoff() const {
  return _reader.getDataField<double>(
      ::capnp::bounded<3>() * ::capnp::ELEMENTS);
}

inline double Stamp::Builder::getFoff() {
  return _builder.getDataField<double>(
      ::capnp::bounded<3>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setFoff(double value) {
  _builder.setDataField<double>(
      ::capnp::bounded<3>() * ::capnp::ELEMENTS, value);
}

inline double Stamp::Reader::getTstart() const {
  return _reader.getDataField<double>(
      ::capnp::bounded<4>() * ::capnp::ELEMENTS);
}

inline double Stamp::Builder::getTstart() {
  return _builder.getDataField<double>(
      ::capnp::bounded<4>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setTstart(double value) {
  _builder.setDataField<double>(
      ::capnp::bounded<4>() * ::capnp::ELEMENTS, value);
}

inline double Stamp::Reader::getTsamp() const {
  return _reader.getDataField<double>(
      ::capnp::bounded<5>() * ::capnp::ELEMENTS);
}

inline double Stamp::Builder::getTsamp() {
  return _builder.getDataField<double>(
      ::capnp::bounded<5>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setTsamp(double value) {
  _builder.setDataField<double>(
      ::capnp::bounded<5>() * ::capnp::ELEMENTS, value);
}

inline  ::int32_t Stamp::Reader::getTelescopeId() const {
  return _reader.getDataField< ::int32_t>(
      ::capnp::bounded<12>() * ::capnp::ELEMENTS);
}

inline  ::int32_t Stamp::Builder::getTelescopeId() {
  return _builder.getDataField< ::int32_t>(
      ::capnp::bounded<12>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setTelescopeId( ::int32_t value) {
  _builder.setDataField< ::int32_t>(
      ::capnp::bounded<12>() * ::capnp::ELEMENTS, value);
}

inline  ::int32_t Stamp::Reader::getNumTimesteps() const {
  return _reader.getDataField< ::int32_t>(
      ::capnp::bounded<13>() * ::capnp::ELEMENTS);
}

inline  ::int32_t Stamp::Builder::getNumTimesteps() {
  return _builder.getDataField< ::int32_t>(
      ::capnp::bounded<13>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setNumTimesteps( ::int32_t value) {
  _builder.setDataField< ::int32_t>(
      ::capnp::bounded<13>() * ::capnp::ELEMENTS, value);
}

inline  ::int32_t Stamp::Reader::getNumChannels() const {
  return _reader.getDataField< ::int32_t>(
      ::capnp::bounded<14>() * ::capnp::ELEMENTS);
}

inline  ::int32_t Stamp::Builder::getNumChannels() {
  return _builder.getDataField< ::int32_t>(
      ::capnp::bounded<14>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setNumChannels( ::int32_t value) {
  _builder.setDataField< ::int32_t>(
      ::capnp::bounded<14>() * ::capnp::ELEMENTS, value);
}

inline  ::int32_t Stamp::Reader::getNumPolarizations() const {
  return _reader.getDataField< ::int32_t>(
      ::capnp::bounded<15>() * ::capnp::ELEMENTS);
}

inline  ::int32_t Stamp::Builder::getNumPolarizations() {
  return _builder.getDataField< ::int32_t>(
      ::capnp::bounded<15>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setNumPolarizations( ::int32_t value) {
  _builder.setDataField< ::int32_t>(
      ::capnp::bounded<15>() * ::capnp::ELEMENTS, value);
}

inline  ::int32_t Stamp::Reader::getNumAntennas() const {
  return _reader.getDataField< ::int32_t>(
      ::capnp::bounded<16>() * ::capnp::ELEMENTS);
}

inline  ::int32_t Stamp::Builder::getNumAntennas() {
  return _builder.getDataField< ::int32_t>(
      ::capnp::bounded<16>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setNumAntennas( ::int32_t value) {
  _builder.setDataField< ::int32_t>(
      ::capnp::bounded<16>() * ::capnp::ELEMENTS, value);
}

inline bool Stamp::Reader::hasData() const {
  return !_reader.getPointerField(
      ::capnp::bounded<1>() * ::capnp::POINTERS).isNull();
}
inline bool Stamp::Builder::hasData() {
  return !_builder.getPointerField(
      ::capnp::bounded<1>() * ::capnp::POINTERS).isNull();
}
inline  ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>::Reader Stamp::Reader::getData() const {
  return ::capnp::_::PointerHelpers< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>>::get(_reader.getPointerField(
      ::capnp::bounded<1>() * ::capnp::POINTERS));
}
inline  ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>::Builder Stamp::Builder::getData() {
  return ::capnp::_::PointerHelpers< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>>::get(_builder.getPointerField(
      ::capnp::bounded<1>() * ::capnp::POINTERS));
}
inline void Stamp::Builder::setData( ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>::Reader value) {
  ::capnp::_::PointerHelpers< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>>::set(_builder.getPointerField(
      ::capnp::bounded<1>() * ::capnp::POINTERS), value);
}
inline void Stamp::Builder::setData(::kj::ArrayPtr<const float> value) {
  ::capnp::_::PointerHelpers< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>>::set(_builder.getPointerField(
      ::capnp::bounded<1>() * ::capnp::POINTERS), value);
}
inline  ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>::Builder Stamp::Builder::initData(unsigned int size) {
  return ::capnp::_::PointerHelpers< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>>::init(_builder.getPointerField(
      ::capnp::bounded<1>() * ::capnp::POINTERS), size);
}
inline void Stamp::Builder::adoptData(
    ::capnp::Orphan< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>>&& value) {
  ::capnp::_::PointerHelpers< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>>::adopt(_builder.getPointerField(
      ::capnp::bounded<1>() * ::capnp::POINTERS), kj::mv(value));
}
inline ::capnp::Orphan< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>> Stamp::Builder::disownData() {
  return ::capnp::_::PointerHelpers< ::capnp::List<float,  ::capnp::Kind::PRIMITIVE>>::disown(_builder.getPointerField(
      ::capnp::bounded<1>() * ::capnp::POINTERS));
}

inline bool Stamp::Reader::hasSeticoreVersion() const {
  return !_reader.getPointerField(
      ::capnp::bounded<2>() * ::capnp::POINTERS).isNull();
}
inline bool Stamp::Builder::hasSeticoreVersion() {
  return !_builder.getPointerField(
      ::capnp::bounded<2>() * ::capnp::POINTERS).isNull();
}
inline  ::capnp::Text::Reader Stamp::Reader::getSeticoreVersion() const {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::get(_reader.getPointerField(
      ::capnp::bounded<2>() * ::capnp::POINTERS));
}
inline  ::capnp::Text::Builder Stamp::Builder::getSeticoreVersion() {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::get(_builder.getPointerField(
      ::capnp::bounded<2>() * ::capnp::POINTERS));
}
inline void Stamp::Builder::setSeticoreVersion( ::capnp::Text::Reader value) {
  ::capnp::_::PointerHelpers< ::capnp::Text>::set(_builder.getPointerField(
      ::capnp::bounded<2>() * ::capnp::POINTERS), value);
}
inline  ::capnp::Text::Builder Stamp::Builder::initSeticoreVersion(unsigned int size) {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::init(_builder.getPointerField(
      ::capnp::bounded<2>() * ::capnp::POINTERS), size);
}
inline void Stamp::Builder::adoptSeticoreVersion(
    ::capnp::Orphan< ::capnp::Text>&& value) {
  ::capnp::_::PointerHelpers< ::capnp::Text>::adopt(_builder.getPointerField(
      ::capnp::bounded<2>() * ::capnp::POINTERS), kj::mv(value));
}
inline ::capnp::Orphan< ::capnp::Text> Stamp::Builder::disownSeticoreVersion() {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::disown(_builder.getPointerField(
      ::capnp::bounded<2>() * ::capnp::POINTERS));
}

inline  ::int32_t Stamp::Reader::getCoarseChannel() const {
  return _reader.getDataField< ::int32_t>(
      ::capnp::bounded<17>() * ::capnp::ELEMENTS);
}

inline  ::int32_t Stamp::Builder::getCoarseChannel() {
  return _builder.getDataField< ::int32_t>(
      ::capnp::bounded<17>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setCoarseChannel( ::int32_t value) {
  _builder.setDataField< ::int32_t>(
      ::capnp::bounded<17>() * ::capnp::ELEMENTS, value);
}

inline  ::int32_t Stamp::Reader::getFftSize() const {
  return _reader.getDataField< ::int32_t>(
      ::capnp::bounded<18>() * ::capnp::ELEMENTS);
}

inline  ::int32_t Stamp::Builder::getFftSize() {
  return _builder.getDataField< ::int32_t>(
      ::capnp::bounded<18>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setFftSize( ::int32_t value) {
  _builder.setDataField< ::int32_t>(
      ::capnp::bounded<18>() * ::capnp::ELEMENTS, value);
}

inline  ::int32_t Stamp::Reader::getStartChannel() const {
  return _reader.getDataField< ::int32_t>(
      ::capnp::bounded<19>() * ::capnp::ELEMENTS);
}

inline  ::int32_t Stamp::Builder::getStartChannel() {
  return _builder.getDataField< ::int32_t>(
      ::capnp::bounded<19>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setStartChannel( ::int32_t value) {
  _builder.setDataField< ::int32_t>(
      ::capnp::bounded<19>() * ::capnp::ELEMENTS, value);
}

inline bool Stamp::Reader::hasSignal() const {
  return !_reader.getPointerField(
      ::capnp::bounded<3>() * ::capnp::POINTERS).isNull();
}
inline bool Stamp::Builder::hasSignal() {
  return !_builder.getPointerField(
      ::capnp::bounded<3>() * ::capnp::POINTERS).isNull();
}
inline  ::Signal::Reader Stamp::Reader::getSignal() const {
  return ::capnp::_::PointerHelpers< ::Signal>::get(_reader.getPointerField(
      ::capnp::bounded<3>() * ::capnp::POINTERS));
}
inline  ::Signal::Builder Stamp::Builder::getSignal() {
  return ::capnp::_::PointerHelpers< ::Signal>::get(_builder.getPointerField(
      ::capnp::bounded<3>() * ::capnp::POINTERS));
}
#if !CAPNP_LITE
inline  ::Signal::Pipeline Stamp::Pipeline::getSignal() {
  return  ::Signal::Pipeline(_typeless.getPointerField(3));
}
#endif  // !CAPNP_LITE
inline void Stamp::Builder::setSignal( ::Signal::Reader value) {
  ::capnp::_::PointerHelpers< ::Signal>::set(_builder.getPointerField(
      ::capnp::bounded<3>() * ::capnp::POINTERS), value);
}
inline  ::Signal::Builder Stamp::Builder::initSignal() {
  return ::capnp::_::PointerHelpers< ::Signal>::init(_builder.getPointerField(
      ::capnp::bounded<3>() * ::capnp::POINTERS));
}
inline void Stamp::Builder::adoptSignal(
    ::capnp::Orphan< ::Signal>&& value) {
  ::capnp::_::PointerHelpers< ::Signal>::adopt(_builder.getPointerField(
      ::capnp::bounded<3>() * ::capnp::POINTERS), kj::mv(value));
}
inline ::capnp::Orphan< ::Signal> Stamp::Builder::disownSignal() {
  return ::capnp::_::PointerHelpers< ::Signal>::disown(_builder.getPointerField(
      ::capnp::bounded<3>() * ::capnp::POINTERS));
}

inline  ::int32_t Stamp::Reader::getSchan() const {
  return _reader.getDataField< ::int32_t>(
      ::capnp::bounded<20>() * ::capnp::ELEMENTS);
}

inline  ::int32_t Stamp::Builder::getSchan() {
  return _builder.getDataField< ::int32_t>(
      ::capnp::bounded<20>() * ::capnp::ELEMENTS);
}
inline void Stamp::Builder::setSchan( ::int32_t value) {
  _builder.setDataField< ::int32_t>(
      ::capnp::bounded<20>() * ::capnp::ELEMENTS, value);
}

inline bool Stamp::Reader::hasObsid() const {
  return !_reader.getPointerField(
      ::capnp::bounded<4>() * ::capnp::POINTERS).isNull();
}
inline bool Stamp::Builder::hasObsid() {
  return !_builder.getPointerField(
      ::capnp::bounded<4>() * ::capnp::POINTERS).isNull();
}
inline  ::capnp::Text::Reader Stamp::Reader::getObsid() const {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::get(_reader.getPointerField(
      ::capnp::bounded<4>() * ::capnp::POINTERS));
}
inline  ::capnp::Text::Builder Stamp::Builder::getObsid() {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::get(_builder.getPointerField(
      ::capnp::bounded<4>() * ::capnp::POINTERS));
}
inline void Stamp::Builder::setObsid( ::capnp::Text::Reader value) {
  ::capnp::_::PointerHelpers< ::capnp::Text>::set(_builder.getPointerField(
      ::capnp::bounded<4>() * ::capnp::POINTERS), value);
}
inline  ::capnp::Text::Builder Stamp::Builder::initObsid(unsigned int size) {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::init(_builder.getPointerField(
      ::capnp::bounded<4>() * ::capnp::POINTERS), size);
}
inline void Stamp::Builder::adoptObsid(
    ::capnp::Orphan< ::capnp::Text>&& value) {
  ::capnp::_::PointerHelpers< ::capnp::Text>::adopt(_builder.getPointerField(
      ::capnp::bounded<4>() * ::capnp::POINTERS), kj::mv(value));
}
inline ::capnp::Orphan< ::capnp::Text> Stamp::Builder::disownObsid() {
  return ::capnp::_::PointerHelpers< ::capnp::Text>::disown(_builder.getPointerField(
      ::capnp::bounded<4>() * ::capnp::POINTERS));
}


CAPNP_END_HEADER


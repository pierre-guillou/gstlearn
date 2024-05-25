#include "Enum/AEnum.hpp"
#include "gstlearn_export.hpp"

namespace hello
{

template <typename T> class GSTLEARN_EXPORT EIterator
{
public:
  using EMap = std::map<int, T*>;

private:
  friend T;
  EIterator() = delete;
  EIterator(EMap* map) : _stditer(map->begin()), _refmap(map) {}

public:
  const T& operator*() const { return (*(_stditer->second)); }

  bool hasNext() const { return (_stditer != _refmap->end()); }
  const T& toNext() { return (*((_stditer++)->second)); }
  const T& toFront()
  {
    _stditer = _refmap->begin();
    return (*(_stditer->second));
  }

  const T& getEnum() const { return (*(_stditer->second)); }

  int getValue() const { return (_stditer->second->getValue()); }

  const String& getKey() const { return (_stditer->second->getKey()); }

  const String& getDescr() const { return (_stditer->second->getDescr()); }

private:
  typename EMap ::iterator _stditer;
  EMap* _refmap;
};

class GSTLEARN_EXPORT EMorpho : public AEnum
{
  using EIterator = EIterator<EMorpho>;
  using EMap = EIterator::EMap;

public:
  EMorpho() : AEnum(*_default) {}
  EMorpho(int value) : AEnum(fromValue(value)) {}
  EMorpho(const String& key) : AEnum(fromKey(key)) {}
  static size_t getSize() { return _map.size(); }
  static EIterator getIterator();
  static void printAll();
  static VectorString getAllKeys();
  static VectorString getAllDescr();
  static bool existsKey(const String& key);
  static bool existsValue(int value)
  {
    return (_map.find(value) != _map.end());
  }

  static const EMorpho& fromKey(const String& key);
  static const EMorpho& fromValue(int value);
  static std ::vector<EMorpho> fromKeys(const VectorString& keys);
  static std ::vector<EMorpho> fromValues(const VectorInt& values);

private:
  EMorpho(const String& key, int value, const String& descr);
  static EMap _map;
  static EIterator _iterator;
  static const EMorpho* _default;

public:
  enum EEMorpho
  {
    E_UNKNOWN = 0,
    E_THRESH = 1,
    E_NEGATION = 2,
    E_EROSION = 3,
    E_DILATION = 4,
    E_OPEN = 5,
    E_CLOSE = 6,
    E_CC = 7,
    E_CCSIZE = 8,
    E_DISTANCE = 9,
    E_ANGLE = 10,
    E_GRADIENT = 11,
  };
  EEMorpho toEnum() const { return static_cast<EEMorpho>(getValue()); }

  static const EMorpho UNKNOWN;
  static const EMorpho THRESH;
  static const EMorpho NEGATION;
  static const EMorpho EROSION;
  static const EMorpho DILATION;
  static const EMorpho OPEN;
  static const EMorpho CLOSE;
  static const EMorpho CC;
  static const EMorpho CCSIZE;
  static const EMorpho DISTANCE;
  static const EMorpho ANGLE;
  static const EMorpho GRADIENT;
};

EMorpho::EMap EMorpho ::_map = EMap();
EMorpho::EIterator EMorpho ::_iterator = EMorpho::EIterator(&EMorpho ::_map);
const EMorpho* EMorpho ::_default = &EMorpho ::UNKNOWN;

EMorpho ::EMorpho(const String& key, int value, const String& descr)
    : AEnum(key, value, descr)
{
  if (_map.find(value) != _map.end())
    throw("Duplicated item");
  _map[value] = this;
}

EMorpho::EIterator EMorpho ::getIterator()
{
  auto it(_iterator);
  it.toFront();
  return it;
}

void EMorpho ::printAll()
{
  auto it(getIterator());
  while (it.hasNext())
    {
      (*it).printEnum();
      it.toNext();
    }
}

VectorString EMorpho ::getAllKeys()
{
  VectorString keys;
  auto it(getIterator());
  while (it.hasNext())
    {
      keys.push_back((*it).getKey());
      it.toNext();
    }
  return keys;
}

VectorString EMorpho ::getAllDescr()
{
  VectorString descr;
  auto it(getIterator());
  while (it.hasNext())
    {
      descr.push_back((*it).getDescr());
      it.toNext();
    }
  return descr;
}

bool EMorpho ::existsKey(const String& key)
{
  auto it = _map.begin();
  while (it != _map.end())
    {
      if (it->second->getKey() == key)
        return true;
      it++;
    }
  return false;
}

const EMorpho& EMorpho ::fromKey(const String& key)
{
  auto it = _map.begin();
  while (it != _map.end())
    {
      if (it->second->getKey() == toUpper(key))
        return (*(it->second));
      it++;
    }
  std ::cout << "Unknown key " << key << " for enum "
             << "EMorpho" << std ::endl;
  return *_default;
}

const EMorpho& EMorpho ::fromValue(int value)
{
  if (existsValue(value))
    return (*(_map[value]));
  std ::cout << "Unknown value " << value << " for enum "
             << "EMorpho" << std ::endl;
  return *_default;
}

std ::vector<EMorpho> EMorpho ::fromKeys(const VectorString& keys)
{
  std ::vector<EMorpho> vec;
  for (auto v : keys)
    vec.push_back(fromKey(v));
  return vec;
}

std ::vector<EMorpho> EMorpho ::fromValues(const VectorInt& values)
{
  std ::vector<EMorpho> vec;
  for (auto v : values)
    vec.push_back(fromValue(v));
  return vec;
}

const EMorpho EMorpho ::UNKNOWN = EMorpho("UNKNOWN", 0, "Idle");
const EMorpho EMorpho ::THRESH =
    EMorpho("THRESH", 1, "Convert the Input Variable into Binary Image");
const EMorpho EMorpho ::NEGATION =
    EMorpho("NEGATION", 2, "Invert of the Binary Image");
const EMorpho EMorpho ::EROSION =
    EMorpho("EROSION", 3, "Erosion on the Binary Image");
const EMorpho EMorpho ::DILATION =
    EMorpho("DILATION", 4, "Dilation on the Binary Image");
const EMorpho EMorpho ::OPEN =
    EMorpho("OPEN", 5, "Opening (erosion then dilation) on the Binary Image");
const EMorpho EMorpho ::CLOSE =
    EMorpho("CLOSE", 6, "Closing (dilation then erosion) on the Binary Image");
const EMorpho EMorpho ::CC =
    EMorpho("CC", 7, "Connected components (cells assigned Rank of CC)");
const EMorpho EMorpho ::CCSIZE =
    EMorpho("CCSIZE", 8, "Connected components (cells assigned Volume of CC)");
const EMorpho EMorpho ::DISTANCE =
    EMorpho("DISTANCE", 9, "Distance to the pore edge");
const EMorpho EMorpho ::ANGLE = EMorpho(
    "ANGLE", 10, "Angle of the tangent to isovalues of a coloured image");
const EMorpho EMorpho ::GRADIENT =
    EMorpho("GRADIENT", 11, "Gradient components");

}; // namespace hello

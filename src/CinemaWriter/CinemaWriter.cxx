#include "CinemaWriter.h"

#include "H5Export.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <vtkAbstractArray.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDirectory.h>
#include <vtkFieldData.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkStringArray.h>

vtkStandardNewMacro(CinemaWriter);

namespace {

using Row = std::map<std::string, std::string>;

std::string JsonQuote(const std::string& value) {
  std::ostringstream out;
  out << '"';
  for (char c : value) {
    if (c == '"' || c == '\\') { out << '\\'; }
    out << c;
  }
  out << '"';
  return out.str();
}

// Converts a field-data array to one scalar or JSON-array CSV value.
std::string FieldArrayValue(vtkAbstractArray* array) {
  if (!array) { return {}; }

  const vtkIdType tuples = array->GetNumberOfTuples();
  const int components = array->GetNumberOfComponents();
  const vtkIdType count = tuples * components;

  if (auto* strings = vtkStringArray::SafeDownCast(array)) {
    if (count == 1) { return strings->GetValue(0); }
    std::ostringstream out;
    out << '[';
    for (vtkIdType i = 0; i < count; ++i) {
      if (i) { out << ','; }
      out << JsonQuote(strings->GetValue(i));
    }
    out << ']';
    return out.str();
  }

  auto* data = vtkDataArray::SafeDownCast(array);
  if (!data) { return {}; }

  if (count == 1) {
    std::ostringstream out;
    out << std::setprecision(17) << data->GetComponent(0, 0);
    return out.str();
  }

  // A single tuple is represented directly as a vector: [x,y,z].
  if (tuples == 1) {
    std::ostringstream out;
    out << '[';
    for (int component = 0; component < components; ++component) {
      if (component) { out << ','; }
      out << std::setprecision(17) << data->GetComponent(0, component);
    }
    out << ']';
    return out.str();
  }

  std::ostringstream out;
  out << '[';
  for (vtkIdType tuple = 0; tuple < tuples; ++tuple) {
    if (tuple) { out << ','; }
    if (components > 1) { out << '['; }
    for (int component = 0; component < components; ++component) {
      if (component) { out << ','; }
      out << std::setprecision(17) << data->GetComponent(tuple, component);
    }
    if (components > 1) { out << ']'; }
  }
  out << ']';
  return out.str();
}

// Collects field data as manifest columns. Array names are used verbatim.
Row CollectFieldAnnotations(vtkDataObject* object) {
  Row row;
  vtkFieldData* fieldData = object ? object->GetFieldData() : nullptr;
  if (!fieldData) { return row; }

  for (int i = 0; i < fieldData->GetNumberOfArrays(); ++i) {
    vtkAbstractArray* array = fieldData->GetAbstractArray(i);
    const char* name = array ? array->GetName() : nullptr;
    if (name && *name) { row[name] = FieldArrayValue(array); }
  }
  return row;
}

// Computes a deterministic 64-bit FNV-1a hash from field-data names, types,
// shapes and values. Only field data contributes to database identity.
std::string FieldDataHash(vtkDataObject* object) {
  vtkFieldData* fieldData = object ? object->GetFieldData() : nullptr;
  std::vector<std::string> records;
  if (fieldData) {
    records.reserve(fieldData->GetNumberOfArrays());
    for (int i = 0; i < fieldData->GetNumberOfArrays(); ++i) {
      vtkAbstractArray* array = fieldData->GetAbstractArray(i);
      if (!array) { continue; }

      std::ostringstream record;
      record << (array->GetName() ? array->GetName() : "") << '|' << array->GetDataType() << '|'
             << array->GetNumberOfTuples() << '|' << array->GetNumberOfComponents() << '|' << FieldArrayValue(array);
      records.push_back(record.str());
    }
  }

  std::sort(records.begin(), records.end());

  std::uint64_t hash = 14695981039346656037ull;
  for (const std::string& record : records) {
    for (unsigned char c : record) {
      hash ^= c;
      hash *= 1099511628211ull;
    }
    hash ^= static_cast<unsigned char>('\n');
    hash *= 1099511628211ull;
  }

  std::ostringstream out;
  out << std::hex << std::setw(16) << std::setfill('0') << hash;
  return out.str();
}

// Recursively returns vtkImageData leaves from a vtkMultiBlockDataSet.
void CollectImages(vtkDataObject* object, std::vector<vtkImageData*>& images) {
  if (auto* image = vtkImageData::SafeDownCast(object)) {
    images.push_back(image);
    return;
  }

  auto* blocks = vtkMultiBlockDataSet::SafeDownCast(object);
  if (!blocks) {
    if (object) { std::cerr << "CinemaWriter: skipping non-image input " << object->GetClassName() << std::endl; }
    return;
  }

  for (unsigned int i = 0; i < blocks->GetNumberOfBlocks(); ++i) { CollectImages(blocks->GetBlock(i), images); }
}

bool WriteImagePNG(vtkImageData* image, const std::string& path, int compressionLevel) {
  if (!image || path.empty()) { return false; }

  vtkPointData* pointData = image->GetPointData();
  if (!pointData) {
    std::cerr << "CinemaWriter: PNG export requires point data." << std::endl;
    return false;
  }

  vtkDataArray* colorArray = nullptr;
  for (int i = 0; i < pointData->GetNumberOfArrays(); ++i) {
    vtkDataArray* array = pointData->GetArray(i);
    const char* name = array ? array->GetName() : nullptr;
    if (!array || !name) { continue; }

    std::string lowerName(name);
    std::transform(lowerName.begin(), lowerName.end(), lowerName.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const int expectedComponents = lowerName == "rgb" ? 3 : (lowerName == "rgba" ? 4 : 0);
    const int dataType = array->GetDataType();
    const bool supportedType = dataType == VTK_UNSIGNED_CHAR || dataType == VTK_UNSIGNED_SHORT;
    if (expectedComponents > 0 && array->GetNumberOfComponents() == expectedComponents && supportedType) {
      colorArray = array;
      break;
    }
  }

  if (!colorArray) {
    std::cerr << "CinemaWriter: PNG export requires an RGB or RGBA point-data array." << std::endl;
    return false;
  }

  vtkNew<vtkImageData> pngImage;
  pngImage->ShallowCopy(image);
  pngImage->GetPointData()->SetScalars(colorArray);

  vtkNew<vtkPNGWriter> writer;
  writer->SetFileName(path.c_str());
  writer->SetCompressionLevel(std::clamp(compressionLevel, 0, 9));
  writer->SetInputData(pngImage);
  writer->Write();

  return writer->GetErrorCode() == 0;
}

std::string CsvEscape(const std::string& value) {
  if (value.find_first_of(",\"\r\n") == std::string::npos) { return value; }
  std::string escaped = "\"";
  for (char c : value) {
    if (c == '"') { escaped += '"'; }
    escaped += c;
  }
  return escaped + '"';
}

// Reads one RFC4180-style CSV record, including quoted values containing
// embedded newlines.
bool ReadCsvRow(std::istream& input, std::vector<std::string>& values) {
  values.clear();
  std::string value;
  bool quoted = false;
  bool readAnything = false;

  for (char c; input.get(c);) {
    readAnything = true;
    if (quoted && c == '"') {
      if (input.peek() == '"') {
        input.get(c);
        value += '"';
      } else {
        quoted = false;
      }
    } else if (c == '"' && value.empty()) {
      quoted = true;
    } else if (c == ',' && !quoted) {
      values.push_back(value);
      value.clear();
    } else if ((c == '\n' || c == '\r') && !quoted) {
      if (c == '\r' && input.peek() == '\n') { input.get(c); }
      values.push_back(value);
      return true;
    } else {
      value += c;
    }
  }

  if (!readAnything) { return false; }
  values.push_back(value);
  return true;
}

void ReadManifest(const std::string& path, std::vector<std::string>& columns, std::vector<Row>& rows) {
  std::ifstream file(path, std::ios::binary);
  if (!file || !ReadCsvRow(file, columns)) { return; }

  std::vector<std::string> values;
  while (ReadCsvRow(file, values)) {
    Row row;
    for (std::size_t i = 0; i < columns.size() && i < values.size(); ++i) { row[columns[i]] = values[i]; }
    rows.push_back(std::move(row));
  }
}

// Rewrites data.csv with the union of old and new columns. FILE and hash are
// kept first; rows sharing hash are replaced.
bool UpdateManifest(const std::string& path, const std::vector<Row>& newRows) {
  const std::string backup = path + ".bak";

  // Recover an interrupted replacement before reading the current manifest.
  std::ifstream current(path);
  const bool hasCurrent = current.good();
  current.close();
  if (!hasCurrent) {
    std::ifstream saved(backup);
    const bool hasBackup = saved.good();
    saved.close();
    if (hasBackup && std::rename(backup.c_str(), path.c_str()) != 0) { return false; }
  }

  std::vector<std::string> columns;
  std::vector<Row> rows;
  ReadManifest(path, columns, rows);

  std::set<std::string> allColumns(columns.begin(), columns.end());
  for (const Row& row : newRows)
    for (const auto& item : row) { allColumns.insert(item.first); }

  columns.clear();
  columns.push_back("FILE");
  columns.push_back("hash");
  allColumns.erase("FILE");
  allColumns.erase("hash");
  columns.insert(columns.end(), allColumns.begin(), allColumns.end());

  for (const Row& incoming : newRows) {
    const auto idIt = incoming.find("hash");
    auto existing = rows.end();
    if (idIt != incoming.end()) {
      existing = std::find_if(rows.begin(), rows.end(), [&](const Row& row) {
        auto it = row.find("hash");
        return it != row.end() && it->second == idIt->second;
      });
    }
    if (existing == rows.end()) {
      rows.push_back(incoming);
    } else {
      *existing = incoming;
    }
  }

  const std::string temp = path + ".tmp";
  std::ofstream file(temp, std::ios::trunc);
  if (!file) { return false; }

  for (std::size_t i = 0; i < columns.size(); ++i) {
    if (i) { file << ','; }
    file << CsvEscape(columns[i]);
  }
  file << '\n';

  for (const Row& row : rows) {
    for (std::size_t i = 0; i < columns.size(); ++i) {
      if (i) { file << ','; }
      auto it = row.find(columns[i]);
      if (it != row.end()) { file << CsvEscape(it->second); }
    }
    file << '\n';
  }
  file.close();
  if (!file) { return false; }

  // Replace the manifest without deleting the last known-good copy first.
  // The backup makes replacement safe on platforms where rename() cannot
  // overwrite an existing file.
  std::remove(backup.c_str());

  std::ifstream existing(path);
  const bool hadExisting = existing.good();
  existing.close();

  if (hadExisting && std::rename(path.c_str(), backup.c_str()) != 0) {
    std::remove(temp.c_str());
    return false;
  }

  if (std::rename(temp.c_str(), path.c_str()) != 0) {
    if (hadExisting) { std::rename(backup.c_str(), path.c_str()); }
    std::remove(temp.c_str());
    return false;
  }

  if (hadExisting) { std::remove(backup.c_str()); }
  return true;
}

} // namespace

CinemaWriter::CinemaWriter() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

CinemaWriter::~CinemaWriter() = default;

void CinemaWriter::ClearProvenance() { this->Provenance.clear(); }

void CinemaWriter::AddProvenanceEntry(const char* key, const char* value) {
  if (!key || !*key) { return; }
  this->Provenance[key] = value ? value : "";
}

int CinemaWriter::FillInputPortInformation(int port, vtkInformation* info) {
  if (port != 0) { return 0; }
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  return 1;
}

int CinemaWriter::FillOutputPortInformation(int port, vtkInformation* info) {
  if (port != 0) { return 0; }
  info->Set(CinemaAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  return 1;
}

// Captures every image leaf as one database/CSV row and passes the input through.
int CinemaWriter::RequestData(vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector) {
  vtkDataObject* input = vtkDataObject::GetData(inputVector[0], 0);
  vtkDataObject* output = vtkDataObject::GetData(outputVector, 0);
  if (!input || !output || this->OutputDirectory.empty()) { return 0; }

  if (!vtkDirectory::MakeDirectory(this->OutputDirectory.c_str())) {
    std::cerr << "CinemaWriter: cannot create output directory " << this->OutputDirectory << std::endl;
    return 0;
  }

  std::vector<vtkImageData*> images;
  CollectImages(input, images);

  std::vector<Row> rows;
  std::set<std::string> fieldDataColumns;
  std::size_t newFiles = 0;
  std::size_t replacedFiles = 0;
  for (vtkImageData* image : images) {
    Row row = CollectFieldAnnotations(image);
    for (const auto& item : row) { fieldDataColumns.insert(item.first); }
    const std::string id = FieldDataHash(image);
    const bool writePNG = this->Format == 1;
    const std::string fileName = id + (writePNG ? ".png" : ".h5");
    const std::string path = this->OutputDirectory + "/" + fileName;
    std::ifstream existingFile(path, std::ios::binary);
    const bool replacingFile = existingFile.good();

    const bool writeSuccess = writePNG ? WriteImagePNG(image, path, this->CompressionLevel) :
                                         WriteImageHDF5(image, path, this->CompressionLevel);
    if (!writeSuccess) { return 0; }
    if (replacingFile) {
      ++replacedFiles;
    } else {
      ++newFiles;
    }

    for (const auto& item : this->Provenance) { row[item.first] = item.second; }
    row["FILE"] = fileName;
    row["hash"] = id;
    rows.push_back(std::move(row));
  }

  if (!UpdateManifest(this->OutputDirectory + "/data.csv", rows)) { return 0; }

  output->ShallowCopy(input);
  std::cout << "# Writer (" << newFiles << " new files, " << replacedFiles << " replaced files, "
            << fieldDataColumns.size() << " field data columns, " << this->Provenance.size() << " provenance columns)"
            << std::endl;
  return 1;
}

// Creates an empty manifest if no manifest exists yet.
int CinemaWriter::CreateDataCSV() const {
  if (this->OutputDirectory.empty()) { return 0; }
  if (!vtkDirectory::MakeDirectory(this->OutputDirectory.c_str())) { return 0; }

  const std::string path = this->OutputDirectory + "/data.csv";
  std::ifstream existing(path);
  if (existing.good()) { return 1; }
  return UpdateManifest(path, {}) ? 1 : 0;
}

// Deletes the complete database contents and recreates the output directory.
int CinemaWriter::DeleteDatabase() const {
  if (this->OutputDirectory.empty()) { return 0; }

  vtkNew<vtkDirectory> directory;
  if (!directory->Open(this->OutputDirectory.c_str())) {
    // A missing directory is already an empty database.
    return vtkDirectory::MakeDirectory(this->OutputDirectory.c_str()) ? 1 : 0;
  }

  if (!vtkDirectory::DeleteDirectory(this->OutputDirectory.c_str())) {
    std::cerr << "CinemaWriter: failed to delete output directory " << this->OutputDirectory << std::endl;
    return 0;
  }

  if (!vtkDirectory::MakeDirectory(this->OutputDirectory.c_str())) {
    std::cerr << "CinemaWriter: failed to recreate output directory " << this->OutputDirectory << std::endl;
    return 0;
  }

  return 1;
}

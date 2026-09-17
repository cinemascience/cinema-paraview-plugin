// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-License-Identifier: BSD-3-Clause
// Cinema database writer for 2D image data.

#pragma once

#include "CinemaWriterModule.h" // for export macro

#include <CinemaAlgorithm.h>
#include <map>
#include <string>

class CINEMAWRITER_EXPORT CinemaWriter : public CinemaAlgorithm {

public:
  static CinemaWriter* New();
  vtkTypeMacro(CinemaWriter, CinemaAlgorithm);

  vtkSetMacro(OutputDirectory, const std::string&);
  vtkGetMacro(OutputDirectory, std::string);

  vtkSetMacro(CompressionLevel, const int);
  vtkGetMacro(CompressionLevel, int);

  vtkSetMacro(Format, const int);
  vtkGetMacro(Format, int);

  int CreateDataCSV() const;
  int DeleteDatabase() const;

  // Populated by vtkSMCinemaWriterProxy immediately before execution.
  // These methods are invoked remotely through the hidden ServerManager
  // Provenance property.
  void ClearProvenance();
  void AddProvenanceEntry(const char* key, const char* value);

protected:
  CinemaWriter();
  ~CinemaWriter();

  int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector) override;
  int FillInputPortInformation(int port, vtkInformation* info) override;
  int FillOutputPortInformation(int port, vtkInformation* info) override;

private:
  CinemaWriter(const CinemaWriter&) = delete;
  void operator=(const CinemaWriter&) = delete;

  std::string OutputDirectory{""};
  int CompressionLevel{5};
  int Format{0};
  std::map<std::string, std::string> Provenance;
};

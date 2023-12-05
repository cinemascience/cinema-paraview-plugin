#pragma once

#include <pcAlgorithm.h>
#include <pcGridLayoutModule.h>

class PCGRIDLAYOUT_EXPORT pcGridLayout : public pcAlgorithm {

private:
  int ColAxis{1};
  int RowAxis{0};

  double ColGap{0};
  double RowGap{0};

  int NumberOfRows{0};

public:
  static pcGridLayout *New();
  vtkTypeMacro(pcGridLayout, pcAlgorithm);

  vtkSetMacro(ColAxis, int);
  vtkGetMacro(ColAxis, int);

  vtkSetMacro(RowAxis, int);
  vtkGetMacro(RowAxis, int);

  vtkSetMacro(ColGap, double);
  vtkGetMacro(ColGap, double);

  vtkSetMacro(RowGap, double);
  vtkGetMacro(RowGap, double);

  vtkSetMacro(NumberOfRows, int);
  vtkGetMacro(NumberOfRows, int);

protected:
  pcGridLayout();
  ~pcGridLayout() override;

  int CopyObject(vtkDataObject *output, vtkDataObject *input);
  int TranslateObject(vtkDataObject *input,
                      const size_t &colAxis,
                      const size_t &rowAxis,
                      const double &dw,
                      const double &dh);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};

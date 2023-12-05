#pragma once

// VTK Module
#include <pcAlgorithmModule.h>

// VTK Includes
#include <vtkAlgorithm.h>
class vtkDataSet;
class vtkInformation;
class vtkInformationIntegerKey;

#include <iostream>

#define pcSetEnumMacro(name, enumType)                    \
  virtual void Set##name(int _arg) {                      \
    vtkDebugMacro(<< this->GetClassName() << " (" << this \
                  << "): setting " #name " to " << _arg); \
    if(this->name != static_cast<enumType>(_arg)) {       \
      this->name = static_cast<enumType>(_arg);           \
    }     \
      this->Modified();                                   \
  }                                                       \
  vtkSetEnumMacro(name, enumType);

class PCALGORITHM_EXPORT pcAlgorithm : public vtkAlgorithm {

public:
  static pcAlgorithm *New();
  vtkTypeMacro(pcAlgorithm, vtkAlgorithm);

  static int inline printErr(const std::string& s) {
    std::cerr<<s<<std::endl;
    return 1;
  }
  static int inline printMsg(const std::string& s) {
    std::cout<<s<<std::endl;
    return 1;
  }

  /**
   * This key can be used during the FillOutputPortInformation() call to
   * specify that an output port should produce the same data type as a
   * certain input port.
   */
  static vtkInformationIntegerKey *SAME_DATA_TYPE_AS_INPUT_PORT();

  /**
   * This method processes a pipeline request such as
   * vtkDemandDrivenPipeline::REQUEST_DATA or
   * vtkDemandDrivenPipeline::REQUEST_INFORMATION.
   *
   * It is not recommended to override this method in order to be conform
   * to the VTK/TTK pipeline model.
   */
  int ProcessRequest(vtkInformation *request,
                     vtkInformationVector **inputVectors,
                     vtkInformationVector *outputVector) override;

  /**
   * Get the output data object for a port on this algorithm.
   */
  vtkDataSet *GetOutput();
  vtkDataSet *GetOutput(int);

  /**
   * Assign a data object as input. Note that this method does not
   * establish a pipeline connection. Use SetInputConnection() to
   * setup a pipeline connection.
   */
  void SetInputData(vtkDataSet *);
  void SetInputData(int, vtkDataSet *);

  /**
   * Assign a data object as input. Note that this method does not
   * establish a pipeline connection. Use AddInputConnection() to
   * setup a pipeline connection.
   */
  void AddInputData(vtkDataSet *);
  void AddInputData(int, vtkDataSet *);

protected:
  pcAlgorithm();
  ~pcAlgorithm() override;

  /**
   * This method is called during the first pipeline pass in
   * ProcessRequest() to create empty output data objects. The data type of
   * the generated outputs is specified in FillOutputPortInformation().
   *
   * In general it should not be necessary to override this method.
   */
  virtual int RequestDataObject(vtkInformation *request,
                                vtkInformationVector **inputVectors,
                                vtkInformationVector *outputVector);

  /**
   * This method is called during the second pipeline pass in
   * ProcessRequest() to provide lightweight information about the outputs
   * without any lengthy computations. For example, the data extend or the
   * number of available time steps.
   *
   * In general, it should only be necessary to override this method to
   * provide information about new vtkImageData output objects, such as
   * their extend, spacing, and origin.
   */
  virtual int RequestInformation(vtkInformation *request,
                       vtkInformationVector **inputVectors,
                       vtkInformationVector *outputVector) {
    return 1;
  }

  /**
   * This method is called during the third pipeline pass in
   * ProcessRequest() to update time.
   *
   * In general it should not be necessary to override this method.
   */
  virtual int
    RequestUpdateTime(vtkInformation *request,
                      vtkInformationVector **inputVectors,
                      vtkInformationVector *outputVector) {
    return 1;
  }

  /**
   * This method is called during the fourth pipeline pass in
   * ProcessRequest() to set time dependent information.
   *
   * In general it should not be necessary to override this method.
   */
  virtual int RequestUpdateTimeDependentInformation(
    vtkInformation *request,
    vtkInformationVector **inputVectors,
    vtkInformationVector *outputVector) {
    return 1;
  }

  /**
   * This method is called during the fifth pipeline pass in
   * ProcessRequest() to specify which portion of its input is needed to
   * create the portion of its output that a downstream filter requested.
   *
   * In general it should not be necessary to override this method unless
   * the filter supports spatial or temporal streaming.
   */
  virtual int
    RequestUpdateExtent(vtkInformation *request,
                        vtkInformationVector **inputVectors,
                        vtkInformationVector *outputVector) {
    return 1;
  }

  /**
   * This method is called during the sixth pipeline pass in
   * ProcessRequest() to specify which outputs will currently not be
   * generated during a RequestData() call.
   *
   * In general it should not be necessary to override this method.
   */
  virtual int
    RequestDataNotGenerated(vtkInformation *request,
                            vtkInformationVector **inputVectors,
                            vtkInformationVector *outputVector) {
    return 1;
  }

  /**
   * This method is called during the seventh pipeline pass in
   * ProcessRequest() to execute an algorithm and update the so far empty
   * output data objects.
   *
   * This method has to be overridden in order to implement the purpose of
   * the filter.
   */
  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVectors,
                          vtkInformationVector *outputVector) {
    return 1;
  }

  /**
   * This method specifies the required input object data types of the
   * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
   * the port information.
   *
   * This method has to be overridden to specify the required input data
   * types.
   */
  int FillInputPortInformation(int port, vtkInformation *info) override {
    return 0;
  }

  /**
   * This method specifies in the port information the data type of the
   * output objects. It is possible to either explicitly specify a type by
   * adding a vtkDataObject::DATA_TYPE_NAME() key, or to pass a type of an
   * input port to an output port by adding the
   * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key.
   *
   * This method has to be overridden to specify the data types of the
   * outputs.
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override {
    return 0;
  }
};

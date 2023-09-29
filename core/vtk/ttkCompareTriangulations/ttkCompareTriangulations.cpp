#include <ttkCompareTriangulations.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <Freudenthal3D.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkCompareTriangulations);

ttkCompareTriangulations::ttkCompareTriangulations() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkCompareTriangulations::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkCompareTriangulations::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkCompareTriangulations::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // this filter just forwards the input to the output
  auto inputDataSet = vtkImageData::GetData(inputVector[0]);
  auto outputDataSet = vtkImageData::GetData(outputVector, 0);
  outputDataSet->ShallowCopy(inputDataSet);

  // compare triangulations
  auto original_triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  this->preconditionTriangulation(original_triangulation); // implemented in base class
  ttkTypeMacroT(
    original_triangulation->getType(),
    this->printTriangulation<T0>( (T0*)original_triangulation->getData() )
  );

  int dim[3];
  inputDataSet->GetDimensions(dim);
  auto lut_triangulation = new ttk::Freudenthal3D({dim[0],dim[1],dim[2]});
  this->preconditionTriangulation(lut_triangulation);
  this->printTriangulation<ttk::Freudenthal3D>( lut_triangulation );

  // return success
  return 1;
}

#include <ttkMacros.h>
#include <ttkMorseSmaleComplex.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSignedCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkMorseSmaleComplex);

ttkMorseSmaleComplex::ttkMorseSmaleComplex() {
  this->setDebugMsgPrefix("MorseSmaleComplex");
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(4);
}

int ttkMorseSmaleComplex::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkMorseSmaleComplex::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 3) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template <typename vtkArrayType, typename vectorType>
void setArray(vtkArrayType &vtkArray, vectorType &vector) {
  ttkUtils::SetVoidArray(vtkArray, vector.data(), vector.size(), 1);
}

template <typename scalarType, typename triangulationType>
int ttkMorseSmaleComplex::dispatch(
  vtkDataArray *const inputScalars,
  vtkDataArray *const inputOffsets,
  vtkUnstructuredGrid *const outputCriticalPoints,
  vtkUnstructuredGrid *const outputSeparatrices1,
  vtkUnstructuredGrid *const outputSeparatrices2,
  const triangulationType &triangulation) {

  const int dimensionality = triangulation.getCellVertexNumber(0) - 1;

  // critical points
  SimplexId criticalPoints_numberOfPoints{};
  std::vector<float> criticalPoints_points;
  std::vector<char> criticalPoints_points_cellDimensions;
  std::vector<scalarType> criticalPoints_points_cellScalars;
  std::vector<SimplexId> criticalPoints_points_cellIds;
  std::vector<char> criticalPoints_points_isOnBoundary;
  std::vector<SimplexId> criticalPoints_points_PLVertexIdentifiers;
  std::vector<SimplexId> criticalPoints_points_manifoldSize;

  // 1-separatrices
  SimplexId separatrices1_numberOfPoints{};
  SimplexId separatrices1_numberOfCells{};

  // 2-separatrices
  SimplexId separatrices2_numberOfPoints{};
  SimplexId separatrices2_numberOfCells{};

  if(ComputeCriticalPoints) {
    this->setOutputCriticalPoints(
      &criticalPoints_numberOfPoints, &criticalPoints_points,
      &criticalPoints_points_cellDimensions, &criticalPoints_points_cellIds,
      &criticalPoints_points_cellScalars, &criticalPoints_points_isOnBoundary,
      &criticalPoints_points_PLVertexIdentifiers,
      &criticalPoints_points_manifoldSize);
  } else {
    this->setOutputCriticalPoints(
      nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
  }

  this->setOutputSeparatrices1(
    &separatrices1_numberOfPoints, &separatrices1_points,
    &separatrices1_points_smoothingMask, &separatrices1_points_cellDimensions,
    &separatrices1_points_cellIds, &separatrices1_numberOfCells,
    &separatrices1_cells_connectivity, &separatrices1_cells_sourceIds,
    &separatrices1_cells_destinationIds, &separatrices1_cells_separatrixIds,
    &separatrices1_cells_separatrixTypes,
    &separatrices1_cells_separatrixFunctionMaxima,
    &separatrices1_cells_separatrixFunctionMinima,
    &separatrices1_cells_separatrixFunctionDiffs,
    &separatrices1_cells_isOnBoundary);

  this->setOutputSeparatrices2(
    &separatrices2_numberOfPoints, &separatrices2_points,
    &separatrices2_numberOfCells, &separatrices2_cells_offsets,
    &separatrices2_cells_connectivity, &separatrices2_cells_sourceIds,
    &separatrices2_cells_separatrixIds, &separatrices2_cells_separatrixTypes,
    &separatrices2_cells_separatrixFunctionMaxima,
    &separatrices2_cells_separatrixFunctionMinima,
    &separatrices2_cells_separatrixFunctionDiffs,
    &separatrices2_cells_isOnBoundary);

  const int ret = this->execute<scalarType, triangulationType>(triangulation);

#ifndef TTK_ENABLE_KAMIKAZE
  if(ret != 0) {
    this->printErr("MorseSmaleComplex.execute() error");
    return -1;
  }
#endif

  // critical points
  {
    vtkNew<vtkPoints> points{};
    vtkNew<vtkSignedCharArray> cellDimensions{};
    vtkNew<ttkSimplexIdTypeArray> cellIds{};
    vtkSmartPointer<vtkDataArray> cellScalars{inputScalars->NewInstance()};
    vtkNew<vtkSignedCharArray> isOnBoundary{};
    vtkNew<ttkSimplexIdTypeArray> PLVertexIdentifiers{};
    vtkNew<ttkSimplexIdTypeArray> manifoldSizeScalars{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!points || !cellDimensions || !cellIds || !cellScalars || !isOnBoundary
       || !PLVertexIdentifiers || !manifoldSizeScalars) {
      this->printErr("Critical points vtkDataArray allocation problem.");
      return -1;
    }
#endif
    points->SetNumberOfPoints(criticalPoints_numberOfPoints);

    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");
    cellDimensions->SetNumberOfTuples(criticalPoints_numberOfPoints);

    cellIds->SetNumberOfComponents(1);
    cellIds->SetName("CellId");
    cellIds->SetNumberOfTuples(criticalPoints_numberOfPoints);

    cellScalars->SetNumberOfComponents(1);
    cellScalars->SetName(inputScalars->GetName());
    cellScalars->SetNumberOfTuples(criticalPoints_numberOfPoints);

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("IsOnBoundary");
    isOnBoundary->SetNumberOfTuples(criticalPoints_numberOfPoints);

    PLVertexIdentifiers->SetNumberOfComponents(1);
    PLVertexIdentifiers->SetName(ttk::VertexScalarFieldName);
    PLVertexIdentifiers->SetNumberOfTuples(criticalPoints_numberOfPoints);

    manifoldSizeScalars->SetNumberOfComponents(1);
    manifoldSizeScalars->SetName("ManifoldSize");
    manifoldSizeScalars->SetNumberOfTuples(criticalPoints_numberOfPoints);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < criticalPoints_numberOfPoints; ++i) {
      points->SetPoint(i, criticalPoints_points[3 * i],
                       criticalPoints_points[3 * i + 1],
                       criticalPoints_points[3 * i + 2]);

      cellDimensions->SetTuple1(i, criticalPoints_points_cellDimensions[i]);
      cellIds->SetTuple1(i, criticalPoints_points_cellIds[i]);
      cellScalars->SetTuple1(i, criticalPoints_points_cellScalars[i]);
      isOnBoundary->SetTuple1(i, criticalPoints_points_isOnBoundary[i]);
      PLVertexIdentifiers->SetTuple1(
        i, criticalPoints_points_PLVertexIdentifiers[i]);

      if(ComputeAscendingSegmentation and ComputeDescendingSegmentation) {
        manifoldSizeScalars->SetTuple1(
          i, criticalPoints_points_manifoldSize[i]);
      } else {
        manifoldSizeScalars->SetTuple1(i, -1);
      }
    }

    outputCriticalPoints->SetPoints(points);

    auto pointData = outputCriticalPoints->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData) {
      this->printErr("outputCriticalPoints has no point data.");
      return -1;
    }
#endif

    pointData->AddArray(cellDimensions);
    pointData->AddArray(cellIds);
    pointData->AddArray(cellScalars);
    pointData->AddArray(isOnBoundary);
    pointData->AddArray(PLVertexIdentifiers);
    pointData->AddArray(manifoldSizeScalars);
  }

  // 1-separatrices
  if(ComputeAscendingSeparatrices1 or ComputeDescendingSeparatrices1
     or ComputeSaddleConnectors) {

    vtkNew<vtkFloatArray> pointsCoords{};
    vtkNew<vtkSignedCharArray> smoothingMask{};
    vtkNew<vtkSignedCharArray> cellDimensions{};
    vtkNew<ttkSimplexIdTypeArray> cellIds{};
    vtkNew<ttkSimplexIdTypeArray> sourceIds{};
    vtkNew<ttkSimplexIdTypeArray> destinationIds{};
    vtkNew<ttkSimplexIdTypeArray> separatrixIds{};
    vtkNew<vtkSignedCharArray> separatrixTypes{};
    vtkNew<vtkDoubleArray> separatrixFunctionMaxima{};
    vtkNew<vtkDoubleArray> separatrixFunctionMinima{};
    vtkNew<vtkDoubleArray> separatrixFunctionDiffs{};
    vtkNew<vtkSignedCharArray> isOnBoundary{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointsCoords || !smoothingMask || !cellDimensions || !cellIds
       || !sourceIds || !destinationIds || !separatrixIds || !separatrixTypes
       || !separatrixFunctionMaxima || !separatrixFunctionMinima
       || !separatrixFunctionDiffs || !isOnBoundary) {
      this->printErr("1-separatrices vtkDataArray allocation problem.");
      return -1;
    }
#endif

    pointsCoords->SetNumberOfComponents(3);
    setArray(pointsCoords, separatrices1_points);

    smoothingMask->SetNumberOfComponents(1);
    smoothingMask->SetName(ttk::MaskScalarFieldName);
    setArray(smoothingMask, separatrices1_points_smoothingMask);

    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");
    setArray(cellDimensions, separatrices1_points_cellDimensions);

    cellIds->SetNumberOfComponents(1);
    cellIds->SetName("CellId");
    setArray(cellIds, separatrices1_points_cellIds);

    sourceIds->SetNumberOfComponents(1);
    sourceIds->SetName("SourceId");
    setArray(sourceIds, separatrices1_cells_sourceIds);

    destinationIds->SetNumberOfComponents(1);
    destinationIds->SetName("DestinationId");
    setArray(destinationIds, separatrices1_cells_destinationIds);

    separatrixIds->SetNumberOfComponents(1);
    separatrixIds->SetName("SeparatrixId");
    setArray(separatrixIds, separatrices1_cells_separatrixIds);

    separatrixTypes->SetNumberOfComponents(1);
    separatrixTypes->SetName("SeparatrixType");
    setArray(separatrixTypes, separatrices1_cells_separatrixTypes);

    separatrixFunctionMaxima->SetNumberOfComponents(1);
    separatrixFunctionMaxima->SetName("SeparatrixFunctionMaximum");
    setArray(
      separatrixFunctionMaxima, separatrices1_cells_separatrixFunctionMaxima);

    separatrixFunctionMinima->SetNumberOfComponents(1);
    separatrixFunctionMinima->SetName("SeparatrixFunctionMinimum");
    setArray(
      separatrixFunctionMinima, separatrices1_cells_separatrixFunctionMinima);

    separatrixFunctionDiffs->SetNumberOfComponents(1);
    separatrixFunctionDiffs->SetName("SeparatrixFunctionDifference");
    setArray(
      separatrixFunctionDiffs, separatrices1_cells_separatrixFunctionDiffs);

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");
    setArray(isOnBoundary, separatrices1_cells_isOnBoundary);

    vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(separatrices1_numberOfCells + 1);
    connectivity->SetNumberOfComponents(1);
    setArray(connectivity, separatrices1_cells_connectivity);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < separatrices1_numberOfCells + 1; ++i) {
      offsets->SetTuple1(i, 2 * i);
    }

    vtkNew<vtkPoints> points{};
    points->SetData(pointsCoords);
    outputSeparatrices1->SetPoints(points);
    vtkNew<vtkCellArray> cells{};
    cells->SetData(offsets, connectivity);
    outputSeparatrices1->SetCells(VTK_LINE, cells);

    auto pointData = outputSeparatrices1->GetPointData();
    auto cellData = outputSeparatrices1->GetCellData();

#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData || !cellData) {
      this->printErr("outputSeparatrices1 has no point or no cell data.");
      return -1;
    }
#endif

    pointData->AddArray(smoothingMask);
    pointData->AddArray(cellDimensions);
    pointData->AddArray(cellIds);

    cellData->AddArray(sourceIds);
    cellData->AddArray(destinationIds);
    cellData->AddArray(separatrixIds);
    cellData->AddArray(separatrixTypes);
    cellData->AddArray(separatrixFunctionMaxima);
    cellData->AddArray(separatrixFunctionMinima);
    cellData->AddArray(separatrixFunctionDiffs);
    cellData->AddArray(isOnBoundary);
  }

  // 2-separatrices
  if(dimensionality == 3
     and (ComputeAscendingSeparatrices2 or ComputeDescendingSeparatrices2)) {

    vtkNew<vtkFloatArray> pointsCoords{};
    vtkNew<ttkSimplexIdTypeArray> sourceIds{};
    vtkNew<ttkSimplexIdTypeArray> separatrixIds{};
    vtkNew<vtkSignedCharArray> separatrixTypes{};
    vtkNew<vtkDoubleArray> separatrixFunctionMaxima{};
    vtkNew<vtkDoubleArray> separatrixFunctionMinima{};
    vtkNew<vtkDoubleArray> separatrixFunctionDiffs{};
    vtkNew<vtkSignedCharArray> isOnBoundary{};

#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointsCoords || !sourceIds || !separatrixIds || !separatrixTypes
       || !separatrixFunctionMaxima || !separatrixFunctionMinima
       || !separatrixFunctionDiffs || !isOnBoundary) {
      this->printErr("2-separatrices vtkDataArray allocation problem.");
      return -1;
    }
#endif

    pointsCoords->SetNumberOfComponents(3);
    setArray(pointsCoords, separatrices2_points);

    sourceIds->SetNumberOfComponents(1);
    sourceIds->SetName("SourceId");
    setArray(sourceIds, separatrices2_cells_sourceIds);

    separatrixIds->SetNumberOfComponents(1);
    separatrixIds->SetName("SeparatrixId");
    setArray(separatrixIds, separatrices2_cells_separatrixIds);

    separatrixTypes->SetNumberOfComponents(1);
    separatrixTypes->SetName("SeparatrixType");
    setArray(separatrixTypes, separatrices2_cells_separatrixTypes);

    separatrixFunctionMaxima->SetNumberOfComponents(1);
    separatrixFunctionMaxima->SetName("SeparatrixFunctionMaximum");
    setArray(
      separatrixFunctionMaxima, separatrices2_cells_separatrixFunctionMaxima);

    separatrixFunctionMinima->SetNumberOfComponents(1);
    separatrixFunctionMinima->SetName("SeparatrixFunctionMinimum");
    setArray(
      separatrixFunctionMinima, separatrices2_cells_separatrixFunctionMinima);

    separatrixFunctionDiffs->SetNumberOfComponents(1);
    separatrixFunctionDiffs->SetName("SeparatrixFunctionDifference");
    setArray(
      separatrixFunctionDiffs, separatrices2_cells_separatrixFunctionDiffs);

    isOnBoundary->SetNumberOfComponents(1);
    isOnBoundary->SetName("NumberOfCriticalPointsOnBoundary");
    setArray(isOnBoundary, separatrices2_cells_isOnBoundary);

    vtkNew<vtkUnsignedCharArray> cellTypes{};
    cellTypes->SetNumberOfComponents(1);
    cellTypes->SetNumberOfTuples(separatrices2_numberOfCells);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < separatrices2_numberOfCells; ++i) {
      if(separatrices2_cells_separatrixTypes[i] == 2) {
        cellTypes->SetTuple1(i, VTK_TRIANGLE);
      } else if(separatrices2_cells_separatrixTypes[i] == 1) {
        cellTypes->SetTuple1(i, VTK_POLYGON);
      }
    }

    vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
    offsets->SetNumberOfComponents(1);
    setArray(offsets, separatrices2_cells_offsets);
    connectivity->SetNumberOfComponents(1);
    setArray(connectivity, separatrices2_cells_connectivity);

    vtkNew<vtkPoints> points{};
    points->SetData(pointsCoords);
    outputSeparatrices2->SetPoints(points);
    vtkNew<vtkCellArray> cells{};
    cells->SetData(offsets, connectivity);
    outputSeparatrices2->SetCells(cellTypes, cells);

    auto cellData = outputSeparatrices2->GetCellData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!cellData) {
      this->printErr("outputSeparatrices2 has no cell data.");
      return -1;
    }
#endif

    cellData->AddArray(sourceIds);
    cellData->AddArray(separatrixIds);
    cellData->AddArray(separatrixTypes);
    cellData->AddArray(separatrixFunctionMaxima);
    cellData->AddArray(separatrixFunctionMinima);
    cellData->AddArray(separatrixFunctionDiffs);
    cellData->AddArray(isOnBoundary);
  }

  return ret;
}

int ttkMorseSmaleComplex::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  const auto input = vtkDataSet::GetData(inputVector[0]);
  auto outputCriticalPoints = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputSeparatrices1 = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto outputSeparatrices2 = vtkUnstructuredGrid::GetData(outputVector, 2);
  auto outputMorseComplexes = vtkDataSet::GetData(outputVector, 3);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    this->printErr("Input pointer is NULL.");
    return -1;
  }
  if(input->GetNumberOfPoints() == 0) {
    this->printErr("Input has no point.");
    return -1;
  }
  if(!outputCriticalPoints or !outputSeparatrices1 or !outputSeparatrices2
     or !outputMorseComplexes) {
    this->printErr("Output pointers are NULL.");
    return -1;
  }
#endif

  const auto triangulation = ttkAlgorithm::GetTriangulation(input);
  if(triangulation == nullptr) {
    this->printErr("Triangulation is null");
    return 0;
  }
  this->preconditionTriangulation(triangulation);

  const auto inputScalars = this->GetInputArrayToProcess(0, inputVector);

#ifndef TTK_ENABLE_KAMIKAZE
  if(inputScalars == nullptr) {
    this->printErr("wrong scalars.");
    return -1;
  }
#endif

  auto inputOffsets = ttkAlgorithm::GetOrderArray(
    input, 0, 1, this->ForceInputOffsetScalarField);

#ifndef TTK_ENABLE_KAMIKAZE
  if(inputOffsets == nullptr) {
    this->printErr("wrong offsets.");
    return -1;
  }
  if(inputOffsets->GetDataType() != VTK_INT
     and inputOffsets->GetDataType() != VTK_ID_TYPE) {
    this->printErr("input offset field type not supported.");
    return -1;
  }
#endif

  this->printMsg("Launching computation on field `"
                 + std::string(inputScalars->GetName()) + "'...");

  // morse complexes
  const SimplexId numberOfVertices = triangulation->getNumberOfVertices();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfVertices) {
    this->printErr("Input has no vertices.");
    return -1;
  }
#endif

  vtkNew<ttkSimplexIdTypeArray> ascendingManifold{};
  vtkNew<ttkSimplexIdTypeArray> descendingManifold{};
  vtkNew<ttkSimplexIdTypeArray> morseSmaleManifold{};
#ifndef TTK_ENABLE_KAMIKAZE
  if(!ascendingManifold || !descendingManifold || !morseSmaleManifold) {
    this->printErr("Manifold vtkDataArray allocation problem.");
    return -1;
  }
#endif
  ascendingManifold->SetNumberOfComponents(1);
  ascendingManifold->SetNumberOfTuples(numberOfVertices);
  ascendingManifold->SetName("AscendingManifold");

  descendingManifold->SetNumberOfComponents(1);
  descendingManifold->SetNumberOfTuples(numberOfVertices);
  descendingManifold->SetName("DescendingManifold");

  morseSmaleManifold->SetNumberOfComponents(1);
  morseSmaleManifold->SetNumberOfTuples(numberOfVertices);
  morseSmaleManifold->SetName("MorseSmaleManifold");

  this->setIterationThreshold(IterationThreshold);

  this->setComputeAscendingSeparatrices1(ComputeAscendingSeparatrices1);
  this->setComputeDescendingSeparatrices1(ComputeDescendingSeparatrices1);
  this->setComputeSaddleConnectors(ComputeSaddleConnectors);

  this->setComputeAscendingSeparatrices2(ComputeAscendingSeparatrices2);
  this->setComputeDescendingSeparatrices2(ComputeDescendingSeparatrices2);

  this->setReturnSaddleConnectors(ReturnSaddleConnectors);
  this->setSaddleConnectorsPersistenceThreshold(
    SaddleConnectorsPersistenceThreshold);

  this->setInputScalarField(ttkUtils::GetVoidPointer(inputScalars));
  this->setInputOffsets(
    static_cast<SimplexId *>(ttkUtils::GetVoidPointer(inputOffsets)));

  void *ascendingManifoldPtr = nullptr;
  void *descendingManifoldPtr = nullptr;
  void *morseSmaleManifoldPtr = nullptr;
  if(ComputeAscendingSegmentation)
    ascendingManifoldPtr = ttkUtils::GetVoidPointer(ascendingManifold);
  if(ComputeDescendingSegmentation)
    descendingManifoldPtr = ttkUtils::GetVoidPointer(descendingManifold);
  if(ComputeAscendingSegmentation and ComputeDescendingSegmentation
     and ComputeFinalSegmentation)
    morseSmaleManifoldPtr = ttkUtils::GetVoidPointer(morseSmaleManifold);

  this->setOutputMorseComplexes(
    ascendingManifoldPtr, descendingManifoldPtr, morseSmaleManifoldPtr);

  int ret{};

  ttkVtkTemplateMacro(
    inputScalars->GetDataType(), triangulation->getType(),
    (ret = dispatch<VTK_TT, TTK_TT>(
       inputScalars, inputOffsets, outputCriticalPoints, outputSeparatrices1,
       outputSeparatrices2, *static_cast<TTK_TT *>(triangulation->getData()))));

  if(ret != 0) {
    return -1;
  }

  outputMorseComplexes->ShallowCopy(input);
  // morse complexes
  if(ComputeAscendingSegmentation or ComputeDescendingSegmentation) {
    vtkPointData *pointData = outputMorseComplexes->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
    if(!pointData) {
      this->printErr("outputMorseComplexes has no point data.");
      return -1;
    }
#endif

    if(ComputeDescendingSegmentation)
      pointData->AddArray(descendingManifold);
    if(ComputeAscendingSegmentation)
      pointData->AddArray(ascendingManifold);
    if(ComputeAscendingSegmentation and ComputeDescendingSegmentation
       and ComputeFinalSegmentation)
      pointData->AddArray(morseSmaleManifold);

    pointData->AddArray(inputOffsets);
  }

  return !ret;
}

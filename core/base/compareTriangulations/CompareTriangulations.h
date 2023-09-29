/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::CompareTriangulations
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %CompareTriangulations class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'CompareTriangulations'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The CompareTriangulations class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class CompareTriangulations : virtual public Debug {

  public:
    CompareTriangulations(){
      this->setDebugMsgPrefix("CompareTriangulations");
    };

    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    template <class triangulationType = ttk::AbstractTriangulation>
    int printTriangulation(const triangulationType *triangulation) const {
      // start global timer
      ttk::Timer globalTimer;

      auto nVertices = triangulation->getNumberOfVertices();
      std::vector<ttk::SimplexId> nNeighbors(nVertices);
      std::vector<std::vector<ttk::SimplexId>> neighbors(nVertices);
      for(ttk::SimplexId v=0; v<nVertices; v++){
        nNeighbors[v] = triangulation->getVertexNeighborNumber(v);
        neighbors[v].resize(nNeighbors[v]);
        for(ttk::SimplexId n=0; n<nNeighbors[v]; n++)
          triangulation->getVertexNeighbor(v,n,neighbors[v][n]);
      }

      this->printMsg(ttk::debug::Separator::L1);
      this->printMsg("nVertices: " + std::to_string(nVertices));
      std::stringstream nNeighborsAsString;
      std::copy(nNeighbors.begin(), nNeighbors.end(), std::ostream_iterator<ttk::SimplexId>(nNeighborsAsString, ","));
      this->printMsg("nNeighbors: [" + nNeighborsAsString.str() + "]");
      for(ttk::SimplexId v=0; v<nVertices; v++){
        std::stringstream neighborsAsString;
        std::copy(neighbors[v].begin(), neighbors[v].end(), std::ostream_iterator<ttk::SimplexId>(neighborsAsString, ","));
        this->printMsg("neighbors["+std::to_string(v)+"]: [" + neighborsAsString.str() + "]");
      }

      this->printMsg(ttk::debug::Separator::L1);
      this->printMsg(
        "Complete", 1, globalTimer.getElapsedTime() // global progress, time
      );
      this->printMsg(ttk::debug::Separator::L1);

      return 1; // return success
    }

  }; // CompareTriangulations class

} // namespace ttk

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
#include <valgrind/callgrind.h>

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
      CALLGRIND_START_INSTRUMENTATION;
      CALLGRIND_TOGGLE_COLLECT;

      // Set which method sets should be compared
      bool testVertexNeigbors = true;
      bool testVertexEdges = true;
      bool testVertexTriangles = true;
      bool testVertexStars = true;

      bool completeOutput = false;
      bool compareImprovement = false;


      ttk::Timer globalTimer;

      /*
       * Fest Vertex Neighbor
       */
      auto nVertices = triangulation->getNumberOfVertices();
      if (completeOutput){
        this->printMsg(ttk::debug::Separator::L1);
        this->printMsg("nVertices: " + std::to_string(nVertices));
      }

      if(testVertexNeigbors){
        //fetch data
        std::vector<ttk::SimplexId> nNeighbors(nVertices);
        std::vector<std::vector<ttk::SimplexId>> neighbors(nVertices);
        for(ttk::SimplexId v=0; v<nVertices; v++){
          nNeighbors[v] = triangulation->getVertexNeighborNumber(v);
          neighbors[v].resize(nNeighbors[v]);
          for(ttk::SimplexId n=0; n<nNeighbors[v]; n++)
            triangulation->getVertexNeighbor(v,n,neighbors[v][n]);
        }

        //print data
        if (completeOutput){
          std::stringstream nNeighborsAsString;
          std::copy(nNeighbors.begin(), nNeighbors.end(), std::ostream_iterator<ttk::SimplexId>(nNeighborsAsString, ","));
          this->printMsg("nNeighbors: [" + nNeighborsAsString.str() + "]");
          for(ttk::SimplexId v=0; v<nVertices; v++){
            std::stringstream neighborsAsString;
            std::copy(neighbors[v].begin(), neighbors[v].end(), std::ostream_iterator<ttk::SimplexId>(neighborsAsString, ","));
            this->printMsg("neighbors["+std::to_string(v)+"]: [" + neighborsAsString.str() + "]");
          }
          //todo timing
        }
      }


      /*
       * Test the edges of vertices
       */
      if(testVertexEdges){
        //fetch data
        std::vector<ttk::SimplexId> nEdges(nVertices);
        std::vector<std::vector<ttk::SimplexId>> edges(nVertices);
        for(ttk::SimplexId v=0; v<nVertices; v++){
          nEdges[v] = compareImprovement? 3 : triangulation->getVertexEdgeNumber(v);
          edges[v].resize(nEdges[v]);
          for(ttk::SimplexId n=0; n<nEdges[v]; n++)
            triangulation->getVertexEdge(v,n,edges[v][n]);
        }

        if (completeOutput){
          //print data
          std::stringstream nEdgesAsString;
          std::copy(nEdges.begin(), nEdges.end(), std::ostream_iterator<ttk::SimplexId>(nEdgesAsString, ","));
          this->printMsg("nEdges: [" + nEdgesAsString.str() + "]");
          for(ttk::SimplexId v=0; v<nVertices; v++){
            std::stringstream edgesAsString;
            std::copy(edges[v].begin(), edges[v].end(), std::ostream_iterator<ttk::SimplexId>(edgesAsString, ","));
            this->printMsg("edges["+std::to_string(v)+"]: [" + edgesAsString.str() + "]");
          }
          //todo timing
        }
      }

      /*
       * Test the triangles of vertices
       */
      if(testVertexTriangles){
        //fetch data
        std::vector<ttk::SimplexId> nTriangles(nVertices);
        std::vector<std::vector<ttk::SimplexId>> triangles(nVertices);
        for(ttk::SimplexId v=0; v<nVertices; v++){
          nTriangles[v] = compareImprovement ? 3 : triangulation->getVertexTriangleNumber(v);
          triangles[v].resize(nTriangles[v]);
          for(ttk::SimplexId n=0; n<nTriangles[v]; n++)
            triangulation->getVertexTriangle(v,n,triangles[v][n]);
        }

        if (completeOutput){
          //print data
          std::stringstream nTrianglesAsString;
          std::copy(nTriangles.begin(), nTriangles.end(), std::ostream_iterator<ttk::SimplexId>(nTrianglesAsString, ","));
          this->printMsg("nTriangles: [" + nTrianglesAsString.str() + "]");
          for(ttk::SimplexId v=0; v<nVertices; v++){
            std::stringstream trianglesAsString;
            std::copy(triangles[v].begin(), triangles[v].end(), std::ostream_iterator<ttk::SimplexId>(trianglesAsString, ","));
            this->printMsg("triangles["+std::to_string(v)+"]: [" + trianglesAsString.str() + "]");
          }
          //todo timing
        }
      }


      /*
       * Test the star of vertices
       */
      if(testVertexStars){
        //fetch data
        std::vector<ttk::SimplexId> nStars(nVertices);
        std::vector<std::vector<ttk::SimplexId>> stars(nVertices);
        for(ttk::SimplexId v=0; v<nVertices; v++){
          nStars[v] = compareImprovement ? 1 : triangulation->getVertexStarNumber(v);
          stars[v].resize(nStars[v]);
          for(ttk::SimplexId n=0; n<nStars[v]; n++)
            triangulation->getVertexStar(v,n,stars[v][n]);
        }

        if (completeOutput){
          //print data
          std::stringstream nStarsAsString;
          std::copy(nStars.begin(), nStars.end(), std::ostream_iterator<ttk::SimplexId>(nStarsAsString, ","));
          this->printMsg("nStars: [" + nStarsAsString.str() + "]");
          for(ttk::SimplexId v=0; v<nVertices; v++){
            std::stringstream starsAsString;
            std::copy(stars[v].begin(), stars[v].end(), std::ostream_iterator<ttk::SimplexId>(starsAsString, ","));
            this->printMsg("stars["+std::to_string(v)+"]: [" + starsAsString.str() + "]");
          }
          //todo timing
        }
      }


      this->printMsg(ttk::debug::Separator::L1);
      auto name = std::string(typeid(triangulation).name());
      std::string cutName = name.substr(9,name.length()-10);
      std::string msg = "Completed " + cutName;
      this->printMsg(msg,1,globalTimer.getElapsedTime());
      this->printMsg(ttk::debug::Separator::L1);


      CALLGRIND_TOGGLE_COLLECT;
      CALLGRIND_STOP_INSTRUMENTATION;

      return 1; // return success
    }

  }; // CompareTriangulations class

} // namespace ttk

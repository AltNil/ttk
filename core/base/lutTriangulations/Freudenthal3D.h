#include <Debug.h>
#include <AbstractTriangulation.h>
#include <iostream>
#include <vector>

namespace ttk {

class Freudenthal3D final : public ttk::AbstractTriangulation {

  public:

    Freudenthal3D(std::array<int,3> extent_) : extent(extent_){

      // precalculate the implicite index offsets for neighbors
      std::vector<std::vector<int>> deltas;
      for(std::vector<std::array<int,3>> caseId : lutNeghborOffset){
        std::vector<int> shifts;
        for(std::array<int,3> coordDeltas: caseId){
          shifts.push_back(coordDeltas[0]+coordDeltas[1]*extent[0]+(extent[0]*extent[1]*coordDeltas[2]));
        }
        lutIndexOffset.push_back(shifts);
      }
    }

    ttk::SimplexId getNumberOfVertices() const final {
      return this->extent[0]*this->extent[1]*this->extent[2];
    };


    /**
     * Returns the number of neighbors for a specified vertex
     * \param vertexId The Id of the vertex for which the number of neighbors should be searched
     * \return Number vor neigbors the given vertex has
    */
    ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId& vertexId)
     const final {
      // return the neigbor by the caseID of the given vertex
      return lutNumNeighbour3d[getCaseID(vertexId)];
    };

    /**
     * get a specified neighbor for a given vertex
     * \param vertexId The id of the vertex from which the neighbbor should be searched for
     * \param localNeighborId The local nubmer of the neighbor to be returned (Should be between 0 and vertexNeighborNumber-1)
     * \param [out] neighborId The output in which the explicit neighborId should be stored in
     * \return 1 if successfull
     */
    inline int getVertexNeighbor(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId
    ) const final {
      // vertexid + offset of the local neighbor dependant of the case the vertex is in
      neighborId = SimplexId(vertexId+lutIndexOffset[getCaseID(vertexId)][localNeighborId]);
      return 1;
    };

    int preconditionVertexNeighbors() final {
      return 1;
    };


  private:
    /**
     * Transpose the ipmlicite index of an vertex into a case identifyer which denotes
     * on which border the vertex is positioned.
     * \param vertexId The Id of the vertex for which the caseId should be searched
     * \return CaseId as an integer between 0 and 26
     */
    inline int getCaseID(const ttk::SimplexId& vertexId) const{
      int ex = extent[0];
      int ey = extent[1];
      int ez = extent[2];

      // Calculate the coordinates of the given vertexId
      auto xyDim = ex*ey;
      int z = vertexId/xyDim;
      int xy = (vertexId-z*xyDim);
      int y = xy/ex;
      int x = xy-y*ex;

      // case Clasification by coordinates
      int caseID = 0;
      
      caseID += (!bool(x)); // +1 if x == 0
      caseID += (!bool(x-(ex-1)))*2; // +2 if x == extent[0]-1
      caseID += (!bool(y))*3; // +3 if y == 0
      caseID += (!bool(y-(ey-1)))*6; // +6 if y == extent[1]-1
      caseID += (!bool(z))*9; // +9 if z == 0
      caseID += (!bool(z-(ez-1)))*18; // +18 if z == extent[2]-1
      if(caseID>=27) return 26; // can not be the case unless Freudenthal3D is used in 2D scenario
      return caseID;
    };

    /*
    int getVertexId(const std::array<int,3>& coordinates)  const{
      auto result =  coordinates[0]+coordinates[1]*extent[0]+(extent[0]*extent[1]*coordinates[2]);
      return result;
    };*/


    const std::array<int,3> extent;
    const int lutNumNeighbour3d[27] = {14,10,10,10,6,8,10,8,6,10,6,8,8,4,7,6,4,4,10,8,6,6,4,4,8,7,4};
    const std::vector<std::vector<std::array<int,3>>> lutNeghborOffset = {
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 1, 1}, {0, 1, 1}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {0, 0, -1}, {0, -1, -1}, {0, -1, 0}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {1, 0, 0}, {1, 0, -1}, {0, 0, -1}},
      {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 1}, {0, 0, -1}, {1, 0, -1}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {0, 0, -1}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {-1, 0, 0}, {-1, 0, 1}, {0, 0, 1}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 0, 1}},
      {{0, -1, 0}, {-1, 0, 0}, {-1, 0, 1}, {0, 0, 1}, {0, -1, -1}, {0, 0, -1}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}},
      {{0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {0, -1, 0}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {1, 0, 0}},
      {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 1}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}},
      {{0, -1, 0}, {-1, 0, 0}, {-1, 0, 1}, {0, 0, 1}, {1, -1, 0}, {1, 0, 0}},
      {{0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 0, 1}},
      {{0, -1, 0}, {-1, 0, 0}, {-1, 0, 1}, {0, 0, 1}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 1, 0}},
      {{0, 0, -1}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {0, -1, -1}, {0, -1, 0}},
      {{0, 0, -1}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 0, -1}, {1, 0, 0}},
      {{0, 0, -1}, {1, 0, -1}, {1, 0, 0}, {0, 1, 0}},
      {{0, 0, -1}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {-1, 0, 0}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}},
      {{0, -1, -1}, {0, 0, -1}, {0, -1, 0}, {-1, 0, 0}}
    };
    std::vector<std::vector<int>> lutIndexOffset;
  };
}
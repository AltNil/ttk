#include <Debug.h>
#include <AbstractTriangulation.h>
#include <iostream>
#include <vector>

namespace ttk {

class Freudenthal3D final : public ttk::AbstractTriangulation {

  public:

    Freudenthal3D(std::array<int,3> extent_) : extent(extent_){
      std::vector<std::vector<int>> deltas;
      int c = 0;
      for(std::vector<std::array<int,3>> caseId : lutNeghborOffset){
        std::vector<int> shifts;
        for(std::array<int,3> coordDeltas: caseId){
          std::array<long,3> longArray = {long(coordDeltas[0]),long(coordDeltas[1]),long(coordDeltas[2])};
          int offset = getVertexId(longArray);
          shifts.push_back(offset);
        }
        deltas.push_back(shifts);
      }
      lutIndexOffset = deltas;
    }

    ttk::SimplexId getNumberOfVertices() const final {
      return this->extent[0]*this->extent[1]*this->extent[2];
    };

    ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId& vertexId)
     const final {
      auto caseID = getCaseID(vertexId);
      return lutNumNeighbour3d[caseID];
      // return 0;
    };

    int getVertexNeighbor(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId
    ) const final {
      auto vid = getCaseID(vertexId);


      /*auto localNeighbourhood = lutNeghborOffset[vid];
      auto pointCoordinates = getCoordinates(vertexId);
      auto offsets = localNeighbourhood[localNeighborId];
      std::array<long,3> newCoords;
      for (int i = 0; i < 3; i++){
        newCoords[i]=pointCoordinates[i]+long(offsets[i]);
      }
      auto newNeighborId = getVertexId(newCoords);*/
      //neighborId = SimplexId(newNeighborId);
      
      auto caseVec = lutIndexOffset[vid];
      neighborId = SimplexId(vertexId+caseVec[localNeighborId]);
      
      //neighborId = SimplexId(vertexId);
      return 1;
    };

    int preconditionVertexNeighbors() final {
      return 1;
    };

    int getCaseID(const ttk::SimplexId& vertexId) const{
      
      std::array<long,3> coords;
      coords = getCoordinates(vertexId);
      
      // case Clasification
      int caseID = 0;
      for (int i = 0; i < 3; i++){
        if (coords[i] ==0){
          caseID += pow(3,i); // 1*i^3
        } else if (coords[i]==extent[i]-1){
          caseID += 2*pow(3,i); // 2*i^3
        }
      }
      return caseID;
    };

    std::array<long,3> getCoordinates(const ttk::SimplexId& vertexId) const{
      auto xyDim = extent[0]*extent[1];
      
      long z = vertexId/xyDim;
      long xy = (vertexId-z*xyDim);
      long y = xy/extent[0];
      long x = xy-y*extent[0];
      std::array<long,3> result = {x,y,z};
      return result;

    };

    int getVertexId(const std::array<long,3> coordinates)  const{
      auto result =  coordinates[0]+coordinates[1]*extent[0]+(extent[0]*extent[1]*coordinates[2]);
      return result;
    };

  private:


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
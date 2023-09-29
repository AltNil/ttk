// #include <Debug.h>
#include <AbstractTriangulation.h>

namespace ttk {

class Freudenthal3D final : public ttk::AbstractTriangulation {

  public:

    Freudenthal3D(std::array<int,3> extent_) : extent(extent_){
    }

    ttk::SimplexId getNumberOfVertices() const final {
      return this->extent[0]*this->extent[1]*this->extent[2];
    };

    ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId& vertexId) const final {
      return 0;
    };

    int getVertexNeighbor(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId
    ) const final {
      return 1;
    };

    int preconditionVertexNeighbors() final {
      return 1;
    };

  private:
    const std::array<int,3> extent;
};

}

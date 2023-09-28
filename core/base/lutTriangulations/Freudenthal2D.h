// #include <Debug.h>
#include <AbstractTriangulation.h>

namespace ttk {

class Freudenthal2D final : public ttk::AbstractTriangulation {

  public:

    Freudenthal3D(std::array<int,2> extent_){
      this->extent = extent_;
    }

    ttk::SimplexId getNumberOfVertices() const final {
      return 0;
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
    std::array<int,2> extent;
};

}

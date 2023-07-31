#ifndef _SVG_H_
#define _SVG_H_

#include "common.h"
#include <fstream>

namespace svg
{
using namespace Eigen;

class SVG {
public:
  SVG(const std::string &filename, const Vector2i &size);
  ~SVG();

  void openAnimGroup();
  void closeAnimGroup(int begin, int end, real_t framerate = 1.f / 12.f);

  struct PathVertexProperty {
    int vertex_index = -1;
    std::string vertex_types = "";
    real_t vertex_depth = -1;
  };

  void writePolyline(const std::vector<Vector2f> &polyline, double stroke_width,
                     const Color &color, bool segments,
                     const std::vector<PathVertexProperty> &vertex_properties =
                         std::vector<PathVertexProperty>());
  void writeDot(Vector2f const &dot, double dot_radius, Color const &color);
  void writeDots(std::vector<Vector2f> const &dots, double dot_radius,
                 std::vector<Color> const &colors);

private:
  void header();
  void footer();

  Vector2i m_size;
  std::ofstream m_file;
};
}

#endif

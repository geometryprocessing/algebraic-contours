#include "svg.h"
#include <string>

namespace svg
{
using namespace std;

SVG::SVG(const std::string &filename, const Vector2i &size) : m_size(size) {

  m_file.open(filename, ofstream::out);
  if (!m_file.is_open()) {
    std::cout << "ERROR: Unable to open SVG output file " << filename
              << std::endl;
  }

  header();
}

SVG::~SVG() {
  footer();
  m_file.close();
}

void SVG::header() {
  m_file << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>"
         << endl
         << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl
         << " \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
         << "<svg viewBox=\"0 0 " << m_size.x() << " " << m_size.y()
         << "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" << endl;
}

void SVG::footer() { m_file << "</svg>" << endl; }

void SVG::openAnimGroup() { m_file << "<g visibility=\"hidden\">"; }

void SVG::closeAnimGroup(int begin, int end, real_t framerate) {
  m_file << "<animate id=\"" << begin
         << "\" attributeName=\"visibility\" attributeType=\"XML\" "
            "begin=\"";
  if (begin > 0) {
    m_file << (begin - 1) << ".end\"";
  } else {
    m_file << "0s;" << end << ".end\"";
  }
  m_file << " dur=\"" << framerate << "s\" to=\"visible\"/>" << endl;
  m_file << "</g>" << endl;
}

void SVG::writePolyline(
    const std::vector<Vector2f> &polyline, double stroke_width,
    const Color &color, bool segments,
    const std::vector<PathVertexProperty> &vertex_properties) {

  m_file << "<g fill=\"none\" stroke=\"rgb(" << int(color[0] * 255) << ","
         << int(color[1] * 255) << "," << int(color[2] * 255)
         << ")\" stroke-width=\"" << stroke_width << "\">";

  m_file << "<path d=\"";
  bool first = true;
  std::string vertex_indices_str = "";
  std::string vertex_types_str = "";
  std::string vertex_depth_str = "";//"z=\"";
  for (size_t i = 0; i < polyline.size(); ++i) {
    const Vector2f &p = polyline.at(i);
    m_file << (first ? "M" : "L") << p.x() << " " << m_size.y() - p.y() << " ";

    if (i < vertex_properties.size() &&
        vertex_properties.at(i).vertex_index >= 0) {
      int idx = vertex_properties.at(i).vertex_index;
      vertex_indices_str +=
          " index" + std::to_string(i) + "=\"" + std::to_string(idx) + "\" ";
    }
    if (i < vertex_properties.size() &&
        !vertex_properties.at(i).vertex_types.empty()) {
      auto type_str = vertex_properties.at(i).vertex_types;
      vertex_types_str +=
          " type_v" + std::to_string(i) + "=\"" + type_str + "\" ";
    }
    if (i < vertex_properties.size() &&
        vertex_properties.at(i).vertex_depth >= 0) {
      auto depth = vertex_properties.at(i).vertex_depth;
      vertex_depth_str += std::to_string(depth);
      if (i + 1 < vertex_properties.size())
        vertex_depth_str += " ";
    }

    if (segments) {
      i++;
      const Vector2f &p = polyline.at(i);

      if (i < vertex_properties.size() &&
          vertex_properties.at(i).vertex_index >= 0) {
        int idx = vertex_properties.at(i).vertex_index;
        vertex_indices_str +=
            " index" + std::to_string(i) + "=\"" + std::to_string(idx) + "\" ";
      }
      if (i < vertex_properties.size() &&
          !vertex_properties.at(i).vertex_types.empty()) {
        auto type_str = vertex_properties.at(i).vertex_types;
        vertex_types_str +=
            " type_v" + std::to_string(i) + "=\"" + type_str + "\" ";
      }
      if (i < vertex_properties.size() &&
          vertex_properties.at(i).vertex_depth >= 0) {
        auto depth = vertex_properties.at(i).vertex_depth;
        vertex_depth_str += std::to_string(depth);
        if (i + 1 < vertex_properties.size())
          vertex_depth_str += " ";
      }

      m_file << "L" << p.x() << " " << m_size.y() - p.y();
      if (i < polyline.size() - 1) {
        m_file << "\"" << vertex_indices_str << vertex_types_str << "/>" << endl
               << "<path d=\"";
        vertex_indices_str = "";
        vertex_types_str = "";
      }
    } else {
      first = false;
    }
  }
  // FIXME vertex_depth_str += "\" ";
  m_file << "\"" << vertex_indices_str << vertex_types_str << vertex_depth_str
         << "/>" << endl;
  m_file << "</g>" << endl;
}

void SVG::writeDot(Vector2f const &dot, double dot_radius, Color const &color) {
  m_file << "<circle cx=\"" << dot.x() << "\" cy=\"" << m_size.y() - dot.y()
         << "\" r=\"" << dot_radius << "\" stroke=\"rgb("
         << int(color[0] * 255) << "," << int(color[1] * 255) << ","
         << int(color[2] * 255) << ")\"" << " fill=\"rgb("
         << int(color[0] * 255) << "," << int(color[1] * 255) << ","
         << int(color[2] * 255) << ")\""
         << "/>" << std::endl;
}

void SVG::writeDots(std::vector<Vector2f> const &dots, double dot_radius,
                    std::vector<Color> const &colors) {
  assert(colors.size() > 0);
  for (size_t i = 0; i < dots.size(); i++) {
    Color c = colors[0];
    if (colors.size() == dots.size())
      c = colors[i];
    writeDot(dots[i], dot_radius, c);
  }
}
}

#include "header.hpp"
#include "serialization.hpp"
#include "error_msg.hpp"
#include <cstdlib>

using namespace std;

/** The number of spaces in indendation */
static const unsigned padding = 2;

// =============================================================================
// Minimal frame reader (input side)

string read_file_to_string(const string& path)
{
  ifstream ifs(path.c_str(), ios::in | ios::binary);
  if(!ifs.good())
    throw error_msg("can not open init-frame file '", path, "'.");
  stringstream ss;
  ss << ifs.rdbuf();
  return ss.str();
}

vector<double> read_frame_field(const string& content, const string& name)
{
  // The field is written as  "name" : { "type" : ..., "value" : [ ... ] }
  // Locate the key, making sure it is really a key (followed by a colon) and
  // not a coincidental substring of a value or a longer name.
  const string key = "\"" + name + "\"";
  size_t pos = string::npos;
  for(size_t from = 0;;)
  {
    const size_t p = content.find(key, from);
    if(p == string::npos) break;
    size_t q = p + key.size();
    while(q < content.size() && (content[q]==' ' || content[q]=='\t')) ++q;
    if(q < content.size() && content[q]==':') { pos = p; break; }
    from = p + key.size();
  }
  if(pos == string::npos)
    throw error_msg("init-frame: field '", name, "' not found in snapshot.");

  // Find this field's "value" array and its bracket range.
  const size_t vpos = content.find("\"value\"", pos);
  if(vpos == string::npos)
    throw error_msg("init-frame: no 'value' for field '", name, "'.");
  const size_t lb = content.find('[', vpos);
  const size_t rb = lb == string::npos ? string::npos : content.find(']', lb);
  if(lb == string::npos || rb == string::npos)
    throw error_msg("init-frame: malformed value array for field '", name, "'.");

  // Parse the doubles between the brackets (the array holds only numbers).
  vector<double> out;
  out.reserve(1u << 18);
  const char* s = content.c_str() + lb + 1;
  const char* const e = content.c_str() + rb;
  char* endp = nullptr;
  while(s < e)
  {
    while(s < e && (*s==',' || *s==' ' || *s=='\n' || *s=='\t' || *s=='\r')) ++s;
    if(s >= e) break;
    const double v = strtod(s, &endp);
    if(endp == s) break; // no numeric progress: stop
    out.push_back(v);
    s = endp;
  }
  return out;
}

oarchive::oarchive(std::ostream& stream_, string id, unsigned version)
  : stream(stream_)
{
  if(!stream.good()) throw bad_stream();
  // write initial brace
  open_group();
  // add description and version
  add("id", id);
  add("version", version);
  // open the data group
  add_key("data");
  open_group();
}

oarchive::~oarchive()
{
  // close the 'data' group
  close_group();
  // write final brace
  if(stream.good()) stream << endl << '}';
}

void oarchive::indent(unsigned n)
{
  level += n;
}

void oarchive::unindent(unsigned n)
{
  level -= n;
}

void oarchive::open_group(const std::string& obrace)
{
  // open brace
  stream << obrace;
  // indent one level
  indent();
  // next element is first
  first = true;
}

void oarchive::close_group(const string& cbrace)
{
  // unindent
  unindent();
  // new line
  stream << endl << string(padding*level, ' ');
  // close brace
  stream << cbrace;
  // not the first in the list
  first = false;
}

void oarchive::add_key(const string& key)
{
  new_line();
  stream << "\"" << key << "\" : ";
}

void oarchive::new_line()
{
  // add comma
  if(!first) stream << ',';
  else first = false;
  // new line
  stream << endl;
  // indent
  stream << string(padding*level, ' ');
}

template<>
void oarchive::add_element<const char*>(const char* const& t)
{
  std::stringstream ss;
  ss << "\"" << t << "\"";
  stream << ss.str();
}

template<>
void oarchive::add_element<std::string>(const std::string& t)
{
  add_element(t.c_str());
}

#ifdef _WIN32
#include <windows.h>
//#include <direct.h>
//#include <io.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <cassert>

using namespace std;


#define PI 3.14159265358979323846

// The outer diameter of the nozzle, which is used for a tolerance in the computation.
#define NOZZLE_DIAMETER 1.5

class vector3
{
public:
  double x;
  double y;
  double z;

  vector3()
  {}

  vector3(double X, double Y, double Z)
    :x(X), y(Y), z(Z)
  {}

  vector3(const vector3 &b)
  {
    this->x = b.x;
    this->y = b.y;
    this->z = b.z;
  }

  vector3 operator - (const vector3 &vv)
  {
    return vector3(x - vv.x, y - vv.y, z - vv.z);
  }

  double& operator [] (int index)
  {
    if (index == 0)
    {
      return x;
    }

    if (index == 1)
    {
      return y;
    }

    if (index == 2)
    {
      return z;
    }

    throw "out of range";
  }

  vector3 operator + (const vector3 &vv)
  {
    return vector3(x + vv.x, y + vv.y, z + vv.z);
  }

  vector3 operator * (const double scale)
  {
    return vector3(x * scale, y * scale, z * scale);
  }

  vector3 operator / (const double scale)
  {
    return vector3(x / scale, y / scale, z / scale);
  }

  double operator * (const vector3 &vv)
  {
    return x * vv.x + y * vv.y + z * vv.z;
  }

  friend ostream& operator <<(ostream& os, const vector3& vv);

  double Length()
  {
    return (sqrt(x * x + y * y + z * z));
  }

  double XYLength()
  {
    return (sqrt(x * x + y * y));
  }

  vector3 MoveDownHalfRadius(double radius)
  {
    return vector3(x, y, z - radius / 2);
    //return vector3(x, y, z);
  }

  vector3 MoveDownRadius(double radius)
  {
    return vector3(x, y, z - radius);
    //return vector3(x, y, z);
  }

  vector3 MoveUp(double radius)
  {
    return vector3(x, y, z + radius);
    //return vector3(x, y, z);
  }

  vector3 MoveDown(double radius)
  {
    return vector3(x, y, z - radius);
    //return vector3(x, y, z);
  }

  double RotateAngle(vector3 &to_vertex)
  {
    Uniform();
    to_vertex.Uniform();
    return asin(x * to_vertex.y - y * to_vertex.x) * 360 / (2 * PI);
  }

  double GetXAngle()
  {
    return GetRotateAngle(y, z);
  }

  double GetYAngle()
  {
    return GetRotateAngle(z, x);
  }

  double GetZAngle()
  {
    return GetRotateAngle(x, y);
  }

  double GetRotateAngle(double input_x, double input_y)
  {
    double temp_value;
    if (input_x != 0)
    {
      temp_value = atan(input_y / input_x) * 180 / PI;
    }
    else
    {
      if (input_y > 0)
      {
        temp_value = 90;
      }
      else
      {
        temp_value = -90;
      }
    }

    if (input_x < 0)
    {
      temp_value += 180;
    }

    return temp_value;
  }

  void Uniform()
  {
    double length = Length();
    x /= length;
    y /= length;
    z /= length;
  }

  void Reverse()
  {
    x = -x;
    y = -y;
    z = -z;
  }
};

ostream& operator<<(ostream& os, const vector3& vv)
{
  os << "<" << vv.x << "," << vv.y << "," << vv.z << "> ";
  return os;
}

class Line
{
public:

  vector3 start_vertex;
  vector3 end_vertex;

  Line(vector3 _start_vertex, vector3 _end_vertex)
    : start_vertex(_start_vertex), end_vertex(_end_vertex)
  {

  }

  vector3 Evaluate(double parameter)
  {
    return start_vertex + (end_vertex - start_vertex) * parameter;
  }

};

class CuttingPlane
{
public:
  vector3 normal;
  double start_distance;

private:
  vector3 uniform_normal;

public:
  CuttingPlane()
  {}

  CuttingPlane(vector3 Normal, double Start_distance)
    : normal(Normal), start_distance(Start_distance)
  {
    uniform_normal = Normal;
    uniform_normal.Uniform();
  }

  friend ostream& operator<<(ostream& os, const CuttingPlane& vv);

  double SignedDistance(vector3 &test_point)
  {
    double distance = uniform_normal * test_point;
    return distance - start_distance;
  }

  bool DeviationFromPlane(vector3 &start_point, vector3 &end_point, double &distance_tolerance)
  {
    double start_distance = SignedDistance(start_point);
    double end_distance = SignedDistance(end_point);

    if (start_distance * end_distance < 0)
    {
      return false;
    }

    if (start_distance * end_distance > 0)
    {
      double min_distance = min(abs(start_distance), abs(end_distance));

      if (min_distance <= distance_tolerance)
      {
        return false;
      }
      else
      {
        return true;
      }
    }

    return false;
  }

  bool OnPlanePositive(vector3 &start_point, vector3 &end_point)
  {
    double start_distance = SignedDistance(start_point);
    double end_distance = SignedDistance(end_point);

    if (start_distance * end_distance > 0)
    {
      if (start_distance > 0)
      {
        return false;
      }
      else
      {
        return true;
      }
    }

    return false;
  }
};

ostream& operator<<(ostream& os, const CuttingPlane& vv)
{
  os << "plane{\n";
  os << vv.normal.x << "*x";
  if (vv.normal.y >= 0)
  {
    os << "+";
  }
  os << vv.normal.y << "*y";
  if (vv.normal.z >= 0)
  {
    os << "+";
  }
  os << vv.normal.z << "*z";

  os << ",";

  os << vv.start_distance << "\n";
  os << "}\n";
  return os;
}

class UniversealCoordinateSystem
{
public:
  vector3 original_point;
  vector3 x_axis;
  vector3 y_axis;

  UniversealCoordinateSystem()
  {
    original_point = vector3(0, 0, 0);
    x_axis = vector3(-1, 0, 0);
    //y_axis = vector3(0, 0, 1);
    y_axis = vector3(0, 1, 0);
  }
};

typedef UniversealCoordinateSystem UCS;

class Circle
{
public:

  UCS _ucs;
  double radius;
  double start_angle;
  double end_angle;

  Circle(double _radius, double _start_angle, double _end_angle)
    : radius(_radius), start_angle(_start_angle), end_angle(_end_angle)
  {

  }

  // parameter is in [0, 1].
  vector3 Evaluate(double parameter)
  {
    parameter = start_angle + (end_angle - start_angle) * parameter;
    return _ucs.original_point + _ucs.x_axis * radius * cos(parameter) + _ucs.y_axis * radius * sin(parameter);
  }

  vector3 EvaluateTangent(double parameter)
  {
    parameter = start_angle + (end_angle - start_angle) * parameter;
    return _ucs.original_point - _ucs.x_axis * radius * sin(parameter) + _ucs.y_axis * radius * cos(parameter);
  }

};

bool GenerateFineGeometryCuttingPlane(vector3 &start_vertex, vector3 &end_vertex, double half_offset_distance, const double &local_layer_height, double last_half_offset_distance, const double &last_local_layer_height, const double &cylinder_height, bool has_cutting_plane, CuttingPlane &cutting_plane, fstream &out)
{
  bool need_to_cut = false;
  double distance_tolerance = NOZZLE_DIAMETER;

  bool not_null = false;

  if (!has_cutting_plane)
  {
    need_to_cut = false;
  }
  else
  {
    if (cutting_plane.DeviationFromPlane(start_vertex, end_vertex, distance_tolerance))
    {
      if (cutting_plane.OnPlanePositive(start_vertex, end_vertex))
      {
        need_to_cut = false;
      }
      else
      {
        return not_null;
      }
    }
    else
    {
      need_to_cut = true;
    }
  }

  vector3 left_start, left_end, right_start, right_end;
  vector3 link_vertex = end_vertex - start_vertex;
  vector3 offset_vertex(link_vertex.y, -link_vertex.x, 0);
  offset_vertex.Uniform();

  left_start = start_vertex + offset_vertex * half_offset_distance;
  left_end = end_vertex + offset_vertex * half_offset_distance;

  offset_vertex.Reverse();
  right_start = start_vertex + offset_vertex * half_offset_distance;
  right_end = end_vertex + offset_vertex * half_offset_distance;

  //vector3 local_x_axis, local_y_axis;
  //local_x_axis = left_end - left_start;
  //local_y_axis = right_start - left_start;
  //local_x_axis.Uniform();
  //local_y_axis.Uniform();
  //vector3 local_origin = left_start.MoveDownRadius(local_layer_height);
  //vector3 corner_link = right_end - local_origin;
  //double local_x_coordinate = corner_link * local_x_axis;
  //double local_y_coordinate = corner_link * local_y_axis;
  double cone_bottom_radius = local_layer_height / 2;
  double cone_top_radius = last_local_layer_height / 2;
  //vector3 box_corner_vertex(local_x_coordinate, local_y_coordinate, local_layer_height);


  double x_length = (left_end - left_start).Length();
  double y_length = (right_start - left_start).Length();
  vector3 rotate_vector = left_end - left_start;
  //   double x_angle = rotate_vector.GetXAngle();
  //   double y_angle = rotate_vector.GetYAngle();
  double z_angle = rotate_vector.GetZAngle();
  vector3 target_center_point = (start_vertex + end_vertex) * 0.5;
  vector3 initial_box_right_corner(x_length / 2, y_length / 2, 0);
  vector3 initial_box_left_corner(-x_length / 2, -y_length / 2, -local_layer_height);

  if (last_local_layer_height == 0 && local_layer_height == 0)
  {
    return not_null;
  }

  if ((start_vertex - end_vertex).Length() == 0)
  {
    return not_null;
  }

  not_null = true;

  if (need_to_cut)
  {
    out << "intersection{\n";
    out << cutting_plane;
  }

  out << "union{\n";

  out << "cone\n\
{ " << left_start.MoveDownRadius(cone_top_radius) << "," << cone_top_radius << "\n" << left_end.MoveDownRadius(cone_bottom_radius) << "," << cone_bottom_radius << "\n\
}\n";

  out << "cone\n\
{ " << right_start.MoveDownRadius(cone_top_radius) << "," << cone_top_radius << "\n" << right_end.MoveDownRadius(cone_bottom_radius) << "," << cone_bottom_radius << "\n\
}\n";

  out << "box\n\
{" << initial_box_left_corner << "," << initial_box_right_corner << "\n\
rotate z*" << z_angle << "\n\
translate " << target_center_point << "\n\
}\n";

  if (last_local_layer_height != 0)
  {
    out << "torus\n\
{" << y_length / 2 << ", " << cone_top_radius << "\n\
rotate x*90 \n\
translate " << (start_vertex).MoveDownRadius(cone_top_radius) << "\n\
}\n";

    vector3 temp = (start_vertex).MoveDownRadius(last_local_layer_height);
    if ((start_vertex - temp).Length() != 0)
    {
      out << "cylinder\n\
{" << (start_vertex).MoveDownRadius(last_local_layer_height) << ", " << (start_vertex) << "," << y_length / 2 << "\n\
}\n";
    }
  }

  if (local_layer_height != 0)
  {
    out << "torus\n\
{" << y_length / 2 << ", " << cone_bottom_radius << "\n\
rotate x*90 \n\
translate " << (end_vertex).MoveDownRadius(cone_bottom_radius) << "\n\
}\n";

    vector3 temp = (end_vertex).MoveDownRadius(local_layer_height);
    if ((end_vertex - temp).Length() != 0)
    {
      out << "cylinder\n\
{" << (end_vertex).MoveDownRadius(local_layer_height) << ", " << (end_vertex) << "," << y_length / 2 << "\n\
}\n";
    }
  }

  if (need_to_cut)
  {
    out << "}\n";
  }

  out << "texture{ GcodeTexture }\n";
  out << "}\n";

  return not_null;
}

void VectorArgumentExtracterStream(const string &argument, int value_start_position, vector3 &vector_value)
{
  try{
    string temp_argument = argument.substr(value_start_position, argument.length() - value_start_position - 1);
    stringstream temp_stream(temp_argument);
    double a, b, c;
    temp_stream >> a >> b >> c;
    vector_value[0] = a;
    vector_value[1] = b;
    vector_value[2] = c;
  }
  catch (exception e)
  {
    return;
  }
}

void ValueArgumentExtracterStream(const string &argument, int value_start_position, double &value)
{
  try{
    string temp_argument = argument.substr(value_start_position, argument.length() - value_start_position);
    stringstream temp_stream(temp_argument);
    double a;
    temp_stream >> a;
    value = a;
  }
  catch (exception e)
  {
    return;
  }
}

void VectorArgumentExtracter(const string &argument, int value_start_position, vector3 &vector_value)
{
  int next_position, current_start = value_start_position, current_count, subfix = 0;
  string temp_sub_string;
  char temp_char;
  double command_value;
  vector3 temp_vector_value = vector_value;

  try{
    for (next_position = value_start_position, current_count = 0; next_position < 100;)
    {
      temp_char = argument[next_position];
      if ((temp_char >= '0' && temp_char <= '9') || temp_char == '.')
      {
        current_count++;
      }
      else if (temp_char == ',')
      {
        temp_sub_string = argument.substr(current_start, current_count);
        temp_vector_value[subfix++] = stod(temp_sub_string);
        current_start = next_position + 1;
        current_count = 0;
      }
      else if (temp_char == ')')
      {
        temp_sub_string = argument.substr(current_start, current_count);
        temp_vector_value[subfix++] = stod(temp_sub_string);
        break;
      }
      next_position++;
    }
  }
  catch (exception e)
  {
    return;
  }

  vector_value = temp_vector_value;
}

void ValueArgumentExtracter(const string &argument, int value_start_position, double &value)
{
  int next_position, current_start = value_start_position, current_count, subfix = 0;
  string temp_sub_string;
  char temp_char;
  double command_value;

  for (next_position = value_start_position, current_count = 0; next_position < 100;)
  {
    temp_char = argument[next_position];
    if ((temp_char >= '0' && temp_char <= '9') || temp_char == '.')
    {
      next_position++;
      current_count++;
    }
    else
    {
      temp_sub_string = argument.substr(current_start, current_count);
      value = stod(temp_sub_string);
      break;
    }
  }
}

double CommandValueExtracter(const string &line_string, int command_position)
{
  int next_position;
  string temp_sub_string;
  char temp_char;
  double command_value;

  for (next_position = command_position + 1; next_position < command_position + 20;)
  {
    temp_char = line_string[next_position];
    if ((temp_char >= '0' && temp_char <= '9') || temp_char == '.' || temp_char == '-')
    {
      next_position++;
    }
    else
    {
      break;
    }
  }

  temp_sub_string = line_string.substr(command_position + 1, next_position - command_position);
  command_value = stod(temp_sub_string);
  return command_value;
}

double FindLocalLayerHeight(std::vector<vector3> *existing_layer, vector3 &test_vertex, bool is_first_layer)
{
  if (is_first_layer)
  {
    return test_vertex.z;
  }

  vector3 temp_vertex, nearest_vertex;
  double minimun_distance = 1000, temp_distance;
  for (int i = 0; i < existing_layer->size(); i++)
  {
    temp_vertex = (*existing_layer)[i];
    temp_distance = (test_vertex - temp_vertex).XYLength();
    if (temp_distance < minimun_distance)
    {
      minimun_distance = temp_distance;
      nearest_vertex = temp_vertex;
    }
  }

  return test_vertex.z - nearest_vertex.z;
}

void Generate_single_file(vector3 camera_position_vector, vector3 look_at, double camera_angle, double nozzle_width, bool has_cutting_plane,
  vector3 filament_color, vector3 cutting_plane_normal, CuttingPlane cutting_plane, double cutting_plane_distance,
  double transparency, double layer_only, double layer_before_and, string in_file_name, bool in_file_name_specified,
  bool look_at_specified, bool angle_specified, bool camera_position_specified, bool layer_only_specified, bool layer_before_and_specified,
  vector3 massive_center, stringstream &out_file_name, int layer_index, bool make_layer_video)
{
  string temp_string, number_string, temp_sub_string;
  char temp_char;
  double temp_number, temp_x, temp_y, temp_z, temp_e, temp_A, temp_B, temp_C;
  double last_e = 0, correct_e;
  double local_layer_height, last_local_layer_height, layer_height_sum, layer_height_power_sum, radius, cylinder_height;
  double current_z;
  double offset_distance, last_offset_distance;
  double aa, bb, cc, delta;
  vector3 last_vertex, current_vertex;
  int e_position, x_position, y_position, z_position, G_position, next_position, temp_end_position, A_position, B_position, C_position;
  int line_count = 0;
  bool has_meet_e = false, first_layer = true;
  double vertices_count = 0;
  double min_x = 10000, min_y = 10000, min_z = 10000, max_x = -10000, max_y = -10000, max_z = -10000;
  bool is_skirt = false;

  double maximum_offset_distance = 0;
  string maximum_offset_command;
  double x_offset = 0;
  double y_offset = 0;

  int layer_count = 0;

  std::vector<vector3> last_layer_vertices, current_layer_vertices;
  std::vector<vector3> *working_layer, *existing_layer, *temp_layer;
  working_layer = &last_layer_vertices;
  existing_layer = &current_layer_vertices;

  if (in_file_name_specified)
  {
    fstream in;
    in.open(in_file_name);

    if (in.is_open())
    {
      // Generate the output file.
      fstream out;
      if (make_layer_video)
        out_file_name << layer_index;
      out.open(out_file_name.str() + ".pov", ofstream::out);

      // Write the header.
      out <<
        "global_settings{\n\
        ambient_light rgb" << filament_color << "*0.1\n\
        assumed_gamma 1\n\
        }\n\
        \n\
        #include \"colors.inc\"\n\
        #include \"textures.inc\"\n\
        \n\
        // A back wall to cast shadows onto\n\
        //plane{ -z, 0\n\
        //pigment{ Gray }\n\
        //finish{ } }\n\
        background{White}\n\
        #declare SpacingX = 20;\n\
        #declare Radius = 5;\n\
        #declare LightX = 15;\n\
        #declare LightY = 40;\n\
        #declare LightZ = -40;\n\
        #declare SRadius = 0;\n\
        #declare SFalloff = 11;\n\
        \n\
        #declare GcodeTexture =\n\
        texture{ pigment{ color rgbf<" << filament_color[0] << ", " << filament_color[1] << ", " << filament_color[2] << ", " << transparency << "> }\n\
        normal{ bumps 0.1 scale 0.15 }\n\
        finish{ phong 0.0 specular 0.0\n\
        ambient rgb" << filament_color << "\n\
        diffuse 0.7\n\
        reflection{ 0.02 metallic 0.5 } }\n\
        } // end of texture\n\
        union{\n";

      // Write the geometry.
      string temp_string, number_string, temp_sub_string;
      char temp_char;
      double temp_number, temp_x, temp_y, temp_z, temp_e, temp_A, temp_B, temp_C;
      double last_e = 0, correct_e;
      double local_layer_height, last_local_layer_height, layer_height_sum, layer_height_power_sum, radius, cylinder_height;
      double current_z;
      double offset_distance, last_offset_distance;
      double aa, bb, cc, delta;
      vector3 last_vertex(0, 0, 0), current_vertex(0, 0, 0);
      int e_position, x_position, y_position, z_position, G_position, next_position, temp_end_position, A_position, B_position, C_position;
      int line_count = 0, layer_count = 0;
      bool has_meet_e = false, first_layer = true;
      double vertices_count = 0;
      double min_x = 10000, min_y = 10000, min_z = 10000, max_x = -10000, max_y = -10000, max_z = -10000;
      bool is_skirt = false;

      double maximum_offset_distance = 0;
      string maximum_offset_command;
      double x_offset = 0;
      double y_offset = 0;

      vertices_count = 0;
      layer_count = 0;
      bool layer_not_null = false;
      while (!in.eof())
      {
        getline(in, temp_string);
        line_count++;

        // Searching for ";(<layer>)".
        if (temp_string.find(";(<layer") == 0)
        {
          layer_count++;
          if (layer_count != 1)
          {
            first_layer = false;
            temp_layer = working_layer;
            working_layer = existing_layer;
            existing_layer = temp_layer;
            working_layer->clear();
            if (layer_not_null)
              out << "texture{ GcodeTexture }\n";
            out << "}\n union {\n";
            layer_not_null = false;
          }
          else
          {
            out << "union {\n";
          }
        }

        e_position = temp_string.find('E');
        x_position = temp_string.find('X');
        y_position = temp_string.find('Y');
        z_position = temp_string.find('Z');
        G_position = temp_string.find('G');
        A_position = temp_string.find('A');
        B_position = temp_string.find('B');
        C_position = temp_string.find('C');
        int semicolon_position = temp_string.find(';');

        if (semicolon_position == 0)
        {
          continue;
        }

        if (G_position != 0)
        {
          // Not a moving command.
          // Continue.
          continue;
        }
        else
        {
          if (temp_string[G_position + 1] == '3')
          {
            is_skirt = true;
          }
          else
          {
            is_skirt = false;
          }
        }

        if (z_position >= 0)
        {
          // Command that changes Z position.
          temp_z = CommandValueExtracter(temp_string, z_position);
          current_z = temp_z;
        }

        if (A_position >= 0)
        {
          // Command that changes Z position.
          temp_A = CommandValueExtracter(temp_string, A_position);
        }

        if (B_position >= 0)
        {
          // Command that changes Z position.
          temp_B = CommandValueExtracter(temp_string, B_position);
        }

        if (C_position >= 0)
        {
          // Command that changes Z position.
          temp_C = CommandValueExtracter(temp_string, C_position);
        }

        if (e_position >= 0)
        {
          // Command that prime filament.
          // Generate cylinder for PovRay to render.

          if (x_position < 0 && y_position < 0 && z_position < 0)
          {
            continue;
          }

          if (!has_meet_e)
          {
            has_meet_e = true;
          }

          temp_e = CommandValueExtracter(temp_string, e_position);

          if (x_position >= 0)
          {
            temp_x = -CommandValueExtracter(temp_string, x_position);
            /*
            Move the model to the original point.
            */
            temp_x += x_offset;
          }
          else
          {
            temp_x = current_vertex.x;
          }

          if (y_position >= 0)
          {
            temp_y = CommandValueExtracter(temp_string, y_position);
            /*
            Move the model to the original point.
            */
            temp_y += y_offset;
          }
          else
          {
            temp_y = current_vertex.y;
          }

          current_vertex.x = temp_x;
          current_vertex.y = temp_y;
          current_vertex.z = current_z;

          cylinder_height = (current_vertex - last_vertex).Length();

          if (cylinder_height == 0)
          {
            continue;
          }

          correct_e = temp_e - last_e;
          last_e = temp_e;

          if (correct_e < 0)
          {
            //last_vertex = current_vertex;
            continue;
          }

          if (correct_e == 0)
          {
            //last_vertex = current_vertex;
            continue;
          }

          local_layer_height = FindLocalLayerHeight(existing_layer, current_vertex, first_layer);
          local_layer_height = max(0.001, local_layer_height);
          offset_distance = (correct_e / (cylinder_height * local_layer_height)) - PI * local_layer_height / 4;

          {
            last_local_layer_height = FindLocalLayerHeight(existing_layer, last_vertex, first_layer);
          }
          last_local_layer_height = max(0.001, last_local_layer_height);
          last_offset_distance = (correct_e / (cylinder_height * last_local_layer_height)) - PI * last_local_layer_height / 4;

          if (local_layer_height == 0 || last_local_layer_height == 0)
          {
            last_vertex = current_vertex;
            continue;
          }

          layer_height_sum = (local_layer_height + last_local_layer_height) / 2;
          layer_height_power_sum = (local_layer_height * local_layer_height + last_local_layer_height * last_local_layer_height) / 4;

          // VOLUME = D * (H + h) / 2 * L + PI * L * ((H / 2)^2 + (h / 2)^2) / 2 + (PI * D * PI * ((H / 2)^2 + (h / 2)^2)) / 4 + PI * (H + h) / 2 * (D / 2)^2.
          // So we have the second order equation according D:
          // D^2 * (PI * (H + h) / 8) + D * ((H + h) / 2 * L + PI^2 * (H^2 + h^2) / 16) + (PI * (H^2 + h^2) * L) / 8 - VOLUME.

          aa = (PI * layer_height_sum / 4.0);
          bb = layer_height_sum * cylinder_height + pow(PI, 2) * layer_height_power_sum / 4;
          cc = (PI * layer_height_power_sum * cylinder_height) / 2 - correct_e;

          delta = pow(bb, 2) - 4 * aa * cc;

          if (delta < 0)
          {
            offset_distance = last_offset_distance = 0;
          }
          else
          {
            last_offset_distance = offset_distance = (-bb + sqrt(delta)) / (2 * aa);
          }

          // We should ignore the volume of the torus in the start and the end. But we keep them in rendering.
          // So the correct offset_distance should be:
          last_offset_distance = offset_distance = (2 * correct_e - layer_height_power_sum * cylinder_height) / (layer_height_sum * cylinder_height * 2);

          // If the width of the path is bigger than the nozzle diameter, force it to be the nozzle diameter.
          if (last_offset_distance > nozzle_width)
          {
            last_offset_distance = offset_distance = nozzle_width;
          }

          if (!is_skirt)
          {
            bool segment_written = false;

            if (make_layer_video)
            {
              if (layer_only_specified)
              {
                if (layer_count == layer_only && layer_count <= layer_index)
                {
                  segment_written = GenerateFineGeometryCuttingPlane(last_vertex, current_vertex, offset_distance / 2, local_layer_height, last_offset_distance / 2, last_local_layer_height, cylinder_height, has_cutting_plane, cutting_plane, out);
                }
              }
              else if (layer_before_and_specified)
              {
                if (layer_count <= layer_before_and && layer_count <= layer_index)
                {
                  segment_written = GenerateFineGeometryCuttingPlane(last_vertex, current_vertex, offset_distance / 2, local_layer_height, last_offset_distance / 2, last_local_layer_height, cylinder_height, has_cutting_plane, cutting_plane, out);
                }
              }
              if (layer_count <= layer_index)
              {
                segment_written = GenerateFineGeometryCuttingPlane(last_vertex, current_vertex, offset_distance / 2, local_layer_height, last_offset_distance / 2, last_local_layer_height, cylinder_height, has_cutting_plane, cutting_plane, out);
              }
            }
            else
            {
              if (layer_only_specified)
              {
                if (layer_count == layer_only)
                {
                  segment_written = GenerateFineGeometryCuttingPlane(last_vertex, current_vertex, offset_distance / 2, local_layer_height, last_offset_distance / 2, last_local_layer_height, cylinder_height, has_cutting_plane, cutting_plane, out);
                }
              }
              else if (layer_before_and_specified)
              {
                if (layer_count <= layer_before_and)
                {
                  segment_written = GenerateFineGeometryCuttingPlane(last_vertex, current_vertex, offset_distance / 2, local_layer_height, last_offset_distance / 2, last_local_layer_height, cylinder_height, has_cutting_plane, cutting_plane, out);
                }
              }
              else
              {
                segment_written = GenerateFineGeometryCuttingPlane(last_vertex, current_vertex, offset_distance / 2, local_layer_height, last_offset_distance / 2, last_local_layer_height, cylinder_height, has_cutting_plane, cutting_plane, out);
              }
            }

            if (segment_written)
            {
              if (!layer_not_null)
              {
                layer_not_null = true;
              }
            }
          }
          last_vertex = current_vertex;

          if (layer_count != 1)
          {
            //mass_center = mass_center + current_vertex;
            vertices_count++;

            if (current_vertex.x < min_x)
            {
              min_x = current_vertex.x;
            }
            if (current_vertex.y < min_y)
            {
              min_y = current_vertex.y;
            }
            if (current_vertex.z < min_z)
            {
              min_z = current_vertex.z;
            }
            if (current_vertex.x > max_x)
            {
              max_x = current_vertex.x;
            }
            if (current_vertex.y > max_y)
            {
              max_y = current_vertex.y;
            }
            if (current_vertex.z > max_z)
            {
              max_z = current_vertex.z;
            }
          }

          // Store the current layer vertices.
          working_layer->push_back(last_vertex);
        }
        else
        {
          if (x_position >= 0)
          {
            if (x_position >= 0)
            {
              temp_x = -CommandValueExtracter(temp_string, x_position);
              /*
              Move the model to the original point.
              */
              temp_x += x_offset;
            }
            else
            {
              temp_x = current_vertex.x;
            }

            if (y_position >= 0)
            {
              temp_y = CommandValueExtracter(temp_string, y_position);
              /*
              Move the model to the original point.
              */
              temp_y += y_offset;
            }
            else
            {
              temp_y = current_vertex.y;
            }

            if (!has_meet_e)
            {
              // First vertex.
              last_vertex.x = temp_x;
              last_vertex.y = temp_y;
              last_vertex.z = current_z;

              continue;
            }
            else
            {
              current_vertex.x = temp_x;
              current_vertex.y = temp_y;
              current_vertex.z = current_z;

              last_vertex = current_vertex;
            }
          }
          else
          {
            if (z_position >= 0)
            {
              // Just moving in Z direction.
              last_vertex.z = current_z;
            }
            else
            {
              continue;
            }
          }
        }
      }

      out << "}\n";
      out << "translate <" << -massive_center.x << ", " << -massive_center.y << ", 0>\n";
      out << "rotate<0,0,-360*(clock+0.00)>\n";
      out << "\n}\n";

      // Recommended look_at.
      if (!look_at_specified)
      {
        look_at = vector3(0, 0, massive_center.z);
      }

      // Recommended camera_angle.
      if (!angle_specified)
      {
        camera_angle = 60;
      }

      // Recommended camera_position.
      if (!camera_position_specified)
      {
        camera_position_vector = look_at;
        double x_edge = max(max_z - min_z, max_x - min_x);
        camera_position_vector.y -= abs((x_edge / 2) / tan(camera_angle / 2) + abs(look_at.y - min_y));
      }

      // Write the light source.
      out <<
"/*light_source {\n\
<0, LightY, LightZ> color White\n\
area_light <8, 0, 0>, <0, 8, 0>, 17, 17\n\
adaptive 0\n\
jitter\n\
\n\
spotlight\n\
point_at <+SpacingX, 0, 0>\n\
tightness 0\n\
radius SRadius\n\
falloff SFalloff\n\
}*/\n\
\n\
//light_source{ <100, 100, 100> color White }\n\
//light_source{ < 500, 1800, -2500> color rgb<1, 1, 1>*0.8 }                // sun\n\
";

      int light_source_count = 0;
      double start_angle, end_angle;
      double bounding_radius = 0;// max(max_z - min_z, (max_x - min_x) / 2);
      double light_source_offset_distance = 30000;  // in millimeters.
      vector3 temp_light_source_position, temp_light_source_tangent;
      double temp_circle_parameter;

      // A quarter of circle.
      start_angle = 0;
      end_angle = 2 * PI;

      Circle circle(bounding_radius + light_source_offset_distance, start_angle, end_angle);
      circle._ucs.original_point[2] = 30000;

      for (int i = 0; i < light_source_count; i++)
      {
        temp_circle_parameter = double(i) / double(light_source_count);
        temp_light_source_position = circle.Evaluate(temp_circle_parameter);
        temp_light_source_tangent = circle.EvaluateTangent(temp_circle_parameter);
        temp_light_source_tangent.Uniform();
        temp_light_source_tangent = temp_light_source_tangent * 3.0;
        out << "light_source{\n" << temp_light_source_position << " color White}\n";
      }

      out << "light_source{ <-3000, 0, 3000> color White }\n\
light_source{ <0, -3000, 3000> color White }\n\
light_source{ <0, 0, 3000> color White }\n\
light_source{ " << camera_position_vector * 10 << "color White*0.5 }\n\
camera{\n\
location " << camera_position_vector << "\n\
right     x*image_width / image_height\n\
angle " << camera_angle << " // direction 2*z\n\
look_at " << look_at << "\n\
sky < 0, 0, 1 >\n\
}\n\
";

      in.close();
      out.close();
    }
  }
}

void GetMassiveCenter(vector3 camera_position_vector, vector3 look_at, double camera_angle, double nozzle_width, bool has_cutting_plane,
  vector3 filament_color, vector3 cutting_plane_normal, CuttingPlane cutting_plane, double cutting_plane_distance,
  double transparency, double layer_only, double layer_before_and, string in_file_name, bool in_file_name_specified,
  bool look_at_specified, bool angle_specified, bool camera_position_specified, bool layer_only_specified, bool layer_before_and_specified,
  vector3 &massive_center/*retval*/, int &layer_count/*retval*/)
{
  string temp_string, number_string, temp_sub_string;
  char temp_char;
  double temp_number, temp_x, temp_y, temp_z, temp_e, temp_A, temp_B, temp_C;
  double last_e = 0, correct_e;
  double local_layer_height, last_local_layer_height, layer_height_sum, layer_height_power_sum, radius, cylinder_height;
  double current_z;
  double offset_distance, last_offset_distance;
  double aa, bb, cc, delta;
  vector3 last_vertex(0, 0, 0), current_vertex(0, 0, 0);
  int e_position, x_position, y_position, z_position, G_position, next_position, temp_end_position, A_position, B_position, C_position;
  int line_count = 0;
  bool has_meet_e = false, first_layer = true;
  double vertices_count = 0;
  double min_x = 10000, min_y = 10000, min_z = 10000, max_x = -10000, max_y = -10000, max_z = -10000;
  bool is_skirt = false;

  double maximum_offset_distance = 0;
  string maximum_offset_command;
  double x_offset = 0;
  double y_offset = 0;

  layer_count = 0;

  std::vector<vector3> last_layer_vertices, current_layer_vertices;
  std::vector<vector3> *working_layer, *existing_layer, *temp_layer;
  working_layer = &last_layer_vertices;
  existing_layer = &current_layer_vertices;

  if (in_file_name_specified)
  {
    ifstream in(in_file_name);
    std::cout << in_file_name << std::endl;
    std::cout << (bool)(in.good()) << std::endl;

    if (in)
    {
      // Compute the messive_center.
      while (in)
      {
        getline(in, temp_string);
        line_count++;

        // Searching for ";(<layer>)".
        if (temp_string.find(";(<layer") == 0)
        {
          layer_count++;
          if (layer_count != 1)
          {
            first_layer = false;
            temp_layer = working_layer;
            working_layer = existing_layer;
            existing_layer = temp_layer;
            working_layer->clear();
            //out << "}\n union {\n";
          }
          else
          {
            //out << "union {\n";
          }
        }

        e_position = temp_string.find('E');
        x_position = temp_string.find('X');
        y_position = temp_string.find('Y');
        z_position = temp_string.find('Z');
        G_position = temp_string.find('G');
        A_position = temp_string.find('A');
        B_position = temp_string.find('B');
        C_position = temp_string.find('C');

        if (G_position != 0)
        {
          // Not a moving command.
          // Continue.
          continue;
        }
        else
        {
          // If you want NOT to output a printing path, change the printing command to G3 ...
          if (temp_string[G_position + 1] == '3')
          {
            is_skirt = true;
          }
          else
          {
            is_skirt = false;
          }
        }

        if (z_position >= 0)
        {
          // Command that changes Z position.
          temp_z = CommandValueExtracter(temp_string, z_position);
          current_z = temp_z;
        }

        if (A_position >= 0)
        {
          // Command that changes Z position.
          temp_A = CommandValueExtracter(temp_string, A_position);
        }

        if (B_position >= 0)
        {
          // Command that changes Z position.
          temp_B = CommandValueExtracter(temp_string, B_position);
        }

        if (C_position >= 0)
        {
          // Command that changes Z position.
          temp_C = CommandValueExtracter(temp_string, C_position);
        }

        if (e_position >= 0)
        {
          // Command that prime filament.
          // Generate cylinder for PovRay to render.

          if (x_position < 0 && y_position < 0 && z_position < 0)
          {
            continue;
          }

          if (!has_meet_e)
          {
            has_meet_e = true;
          }

          temp_e = CommandValueExtracter(temp_string, e_position);

          if (x_position >= 0)
          {
            temp_x = -CommandValueExtracter(temp_string, x_position);
            /*
            Move the model to the original point.
            */
            temp_x += x_offset;
          }
          else
          {
            temp_x = current_vertex.x;
          }

          if (y_position >= 0)
          {
            temp_y = CommandValueExtracter(temp_string, y_position);
            /*
            Move the model to the original point.
            */
            temp_y += y_offset;
          }
          else
          {
            temp_y = current_vertex.y;
          }

          current_vertex.x = temp_x;
          current_vertex.y = temp_y;
          current_vertex.z = current_z;

          cylinder_height = (current_vertex - last_vertex).Length();

          if (cylinder_height == 0)
          {
            continue;
          }

          correct_e = temp_e - last_e;
          last_e = temp_e;

          if (correct_e < 0)
          {
            //last_vertex = current_vertex;
            continue;
          }

          if (correct_e == 0)
          {
            //last_vertex = current_vertex;
            continue;
          }

          //correct_e = PI * 2.25 * correct_e;

          local_layer_height = FindLocalLayerHeight(existing_layer, current_vertex, first_layer);
          local_layer_height = max(0.001, local_layer_height);
          offset_distance = (correct_e / (cylinder_height * local_layer_height)) - PI * local_layer_height / 4;

          last_local_layer_height = FindLocalLayerHeight(existing_layer, last_vertex, first_layer);
          last_local_layer_height = max(0.001, last_local_layer_height);
          last_offset_distance = (correct_e / (cylinder_height * last_local_layer_height)) - PI * last_local_layer_height / 4;

          if (local_layer_height == 0 || last_local_layer_height == 0)
          {
            last_vertex = current_vertex;
            continue;
          }

          layer_height_sum = (local_layer_height + last_local_layer_height) / 2;
          layer_height_power_sum = (local_layer_height * local_layer_height + last_local_layer_height * last_local_layer_height) / 4;

          // VOLUME = D * (H + h) / 2 * L + PI * L * ((H / 2)^2 + (h / 2)^2) / 2 + (PI * D * PI * ((H / 2)^2 + (h / 2)^2)) / 4 + PI * (H + h) / 2 * (D / 2)^2.
          // So we have the second order equation according D:
          // D^2 * (PI * (H + h) / 8) + D * ((H + h) / 2 * L + PI^2 * (H^2 + h^2) / 16) + (PI * (H^2 + h^2) * L) / 8 - VOLUME.

          aa = (PI * layer_height_sum / 4.0);
          bb = layer_height_sum * cylinder_height + pow(PI, 2) * layer_height_power_sum / 4;
          cc = (PI * layer_height_power_sum * cylinder_height) / 2 - correct_e;

          delta = pow(bb, 2) - 4 * aa * cc;

          if (delta < 0)
          {
            offset_distance = last_offset_distance = 0;
          }
          else
          {
            last_offset_distance = offset_distance = (-bb + sqrt(delta)) / (2 * aa);
          }

          // We should ignore the volume of the torus in the start and the end. But we keep them in rendering.
          // So the correct offset_distance should be:
          last_offset_distance = offset_distance = (2 * correct_e - layer_height_power_sum * cylinder_height) / (layer_height_sum * cylinder_height * 2);

          // If the width of the path is bigger than the nozzle diameter, force it to be the nozzle diameter.
          if (last_offset_distance > nozzle_width)
          {
            last_offset_distance = offset_distance = nozzle_width;
          }

          last_vertex = current_vertex;

          if (layer_count != 1)
          {
            massive_center = massive_center + current_vertex;
            vertices_count++;

            if (current_vertex.x < min_x)
            {
              min_x = current_vertex.x;
            }
            if (current_vertex.y < min_y)
            {
              min_y = current_vertex.y;
            }
            if (current_vertex.z < min_z)
            {
              min_z = current_vertex.z;
            }
            if (current_vertex.x > max_x)
            {
              max_x = current_vertex.x;
            }
            if (current_vertex.y > max_y)
            {
              max_y = current_vertex.y;
            }
            if (current_vertex.z > max_z)
            {
              max_z = current_vertex.z;
            }
          }
          // Store the current layer vertices.
          working_layer->push_back(last_vertex);
        }
        else
        {
          if (x_position >= 0)
          {
            if (x_position >= 0)
            {
              temp_x = -CommandValueExtracter(temp_string, x_position);
              /*
              Move the model to the original point.
              */
              temp_x += x_offset;
            }
            else
            {
              temp_x = current_vertex.x;
            }

            if (y_position >= 0)
            {
              temp_y = CommandValueExtracter(temp_string, y_position);
              /*
              Move the model to the original point.
              */
              temp_y += y_offset;
            }
            else
            {
              temp_y = current_vertex.y;
            }

            if (!has_meet_e)
            {
              // First vertex.
              last_vertex.x = temp_x;
              last_vertex.y = temp_y;
              last_vertex.z = current_z;

              continue;
            }
            else
            {
              current_vertex.x = temp_x;
              current_vertex.y = temp_y;
              current_vertex.z = current_z;

              last_vertex = current_vertex;
            }
          }
          else
          {
            if (z_position >= 0)
            {
              // Just moving in Z direction.
              last_vertex.z = current_z;
            }
            else
            {
              continue;
            }
          }
        }
      }
      massive_center = massive_center / vertices_count;
      in.close();
    }
  }
}

int main(int argc, char* argv[])
{
  vector3 camera_position_vector(-100, -250, 220);
  vector3 look_at(-100, 100, 30);
  double camera_angle = 60;
  double nozzle_width = 0.4;
  bool has_cutting_plane = false;
  vector3 filament_color = vector3(0.3, 0.8, 1.0) * 0.59;
  vector3 cutting_plane_normal = vector3(0.3, 0.8, 1.0) * 0.59;
  CuttingPlane cutting_plane;
  double cutting_plane_distance = 0.0;
  double transparency = 0.0;
  double layer_only = -1;
  double layer_before_and = -1;
  string in_file_name;
  bool in_file_name_specified = false;
  bool look_at_specified = false;
  bool angle_specified = false;
  bool camera_position_specified = false;
  bool layer_only_specified = false;
  bool layer_before_and_specified = false;
  bool make_layer_video = false;

  // Read configurations file.
  cout << "Reading configurations file...\n";
  fstream configure_file;
  configure_file.open("configurations.txt");
  if (!configure_file)
  {
    cout << "Configurations file not found.\n" << "Exit.\n";
    exit(0);
  }
  else
  {
    string temp_string;
    cout << "Parameters derived:\n";

    while (!configure_file.eof())
    {
      getline(configure_file, temp_string);

      if (temp_string[0] == '#' && temp_string[1] == '#')
        continue;

      int camera_location_conf_pos = temp_string.find("camera_location=v(");
      int camera_lookat_conf_pos = temp_string.find("camera_lookat=v(");
      int filament_color_conf_pos = temp_string.find("filament_color=rgb(");
      int view_angle_conf_pos = temp_string.find("view_angle=");
      int nozzle_conf_pos = temp_string.find("nozzle=");
      int file_conf_pos = temp_string.find("file=");
      int transparency_conf_pos = temp_string.find("transparency=");
      int layer_only_pos = temp_string.find("layer_only=");
      int layer_before_and_pos = temp_string.find("layer_before_and=");
      int cutting_plane_normal_pos = temp_string.find("cutting_plane_normal=");
      int cutting_plane_distance_pos = temp_string.find("cutting_plane_distance=");
      int make_layer_video_pos = temp_string.find("make_layer_video");

      //////////////////////////////////////////////////////////////////////////
      // Boolean arguments.
      if (make_layer_video_pos >= 0)
      {
        make_layer_video = true;
        cout << "make_layer_video=yes\n";
      }

      //////////////////////////////////////////////////////////////////////////
      // Vector arguments.
      if (camera_location_conf_pos >= 0)
      {
        VectorArgumentExtracterStream(temp_string, camera_location_conf_pos + string("camera_location=v(").length(), camera_position_vector);
        camera_position_specified = true;
        cout << "camera_location=" << camera_position_vector << "\n";
      }

      if (camera_lookat_conf_pos >= 0)
      {
        VectorArgumentExtracterStream(temp_string, camera_lookat_conf_pos + string("camera_lookat=v(").length(), look_at);
        look_at_specified = true;
        cout << "camera_lookat=" << look_at << "\n";
      }

      if (filament_color_conf_pos >= 0)
      {
        VectorArgumentExtracterStream(temp_string, filament_color_conf_pos + string("filament_color=rgb(").length(), filament_color);
        cout << "filament_color=" << filament_color << "\n";
      }

      if (cutting_plane_normal_pos >= 0)
      {
        VectorArgumentExtracterStream(temp_string, cutting_plane_normal_pos + string("cutting_plane_normal=v(").length(), cutting_plane_normal);
        has_cutting_plane = true;
        cout << "cutting_plane_normal=" << cutting_plane_normal << "\n";
      }

      //////////////////////////////////////////////////////////////////////////
      // Value arguments.
      if (view_angle_conf_pos >= 0)
      {
        ValueArgumentExtracterStream(temp_string, view_angle_conf_pos + string("view_angle=").length(), camera_angle);
        angle_specified = true;
        cout << "view_angle=" << camera_angle << "\n";
      }

      if (nozzle_conf_pos >= 0)
      {
        ValueArgumentExtracterStream(temp_string, nozzle_conf_pos + string("nozzle_diameter=").length(), nozzle_width);
        cout << "nozzle_diameter=" << nozzle_width << "\n";
      }

      if (transparency_conf_pos >= 0)
      {
        ValueArgumentExtracterStream(temp_string, transparency_conf_pos + string("transparency=").length(), transparency);
        cout << "transparency=" << transparency << "\n";
      }

      if (layer_only_pos >= 0)
      {
        ValueArgumentExtracterStream(temp_string, layer_only_pos + string("layer_only=").length(), layer_only);
        if (layer_only >= 0)
        {
          layer_only_specified = true;
          cout << "layer_only=" << layer_only << "\n";
        }
      }

      if (layer_before_and_pos >= 0)
      {
        ValueArgumentExtracterStream(temp_string, layer_before_and_pos + string("layer_before_and=").length(), layer_before_and);
        if (layer_before_and >= 0)
        {
          layer_before_and_specified = true;
          cout << "layer_before_and=" << layer_before_and << "\n";
        }
      }

      if (cutting_plane_distance_pos >= 0)
      {
        ValueArgumentExtracterStream(temp_string, cutting_plane_distance_pos + string("cutting_plane_distance=").length(), cutting_plane_distance);
        cout << "cutting_plane_distance=" << cutting_plane_distance << "\n";
      }

      //////////////////////////////////////////////////////////////////////////
      // String arguments.
      if (file_conf_pos >= 0)
      {
        in_file_name = temp_string.substr(file_conf_pos + string("file=").length(), temp_string.size() - file_conf_pos - string("file=").length());
        in_file_name_specified = true;
        istringstream ss(in_file_name);
        ss >> in_file_name;
      }
    }

    if (has_cutting_plane)
    {
      cutting_plane = CuttingPlane(cutting_plane_normal, cutting_plane_distance);
    }
  }

  // Compute the massive center.
  vector3 massive_center(0, 0, 0);
  int layer_count = 0;
  GetMassiveCenter(camera_position_vector, look_at, camera_angle, nozzle_width, has_cutting_plane,
    filament_color, cutting_plane_normal, cutting_plane, cutting_plane_distance,
    transparency, layer_only, layer_before_and, in_file_name, in_file_name_specified,
    look_at_specified, angle_specified, camera_position_specified, layer_only_specified, layer_before_and_specified,
    massive_center, layer_count);

  cout << "Layer count: " <<layer_count << "\n";
  cout << "Gcode Massive center: " << massive_center << "\n";
  cout << "Move the hole print to the original point.\n";

  // Write to the POV-Ray file.
  cout << "Writing POV-Ray file...\n";
  if (make_layer_video)
  {
    stringstream temp;
    temp << in_file_name.substr(0, in_file_name.find(".gcode"));
    // _mkdir(temp.str().c_str());

    for (int layer_index = 0; layer_index < layer_count; layer_index++)
    {
      stringstream out_file_name;
      out_file_name << in_file_name.substr(0, in_file_name.find(".gcode")) << "\\";

      Generate_single_file(camera_position_vector, look_at, camera_angle, nozzle_width, has_cutting_plane,
        filament_color, cutting_plane_normal, cutting_plane, cutting_plane_distance,
        transparency, layer_only, layer_before_and, in_file_name, in_file_name_specified,
        look_at_specified, angle_specified, camera_position_specified, layer_only_specified, layer_before_and_specified,
        massive_center, out_file_name, layer_index, make_layer_video);
    }
  }
  else
  {
    stringstream out_file_name;
    out_file_name << in_file_name.substr(0, in_file_name.find(".gcode"));
#ifdef _WIN32
    SYSTEMTIME sys;
    GetLocalTime(&sys);
    out_file_name << '_' << sys.wMonth << '_' << sys.wDay << '_' << sys.wHour << '_' << sys.wMinute;
#endif

    Generate_single_file(camera_position_vector, look_at, camera_angle, nozzle_width, has_cutting_plane,
      filament_color, cutting_plane_normal, cutting_plane, cutting_plane_distance,
      transparency, layer_only, layer_before_and, in_file_name, in_file_name_specified,
      look_at_specified, angle_specified, camera_position_specified, layer_only_specified, layer_before_and_specified,
      massive_center, out_file_name, 0, make_layer_video);
  }

  return 0;
}

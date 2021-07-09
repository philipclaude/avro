#version 410

layout (triangles, equal_spacing, ccw) in;

out vec3 v_Position;
out vec2 v_ParameterTess;
out vec3 v_Normal;

void main() {

  float s = gl_TessCoord.x;
  float t = gl_TessCoord.y;

  vec3 p0 = gl_in[0].gl_Position.xyz;
  vec3 p1 = gl_in[1].gl_Position.xyz;
  vec3 p2 = gl_in[2].gl_Position.xyz;

  // TODO, use basis functions to calculate normal
  v_Normal = normalize(cross( p1 - p0 , p2 - p0 ));

  #if GEOMETRY_ORDER == 1
  v_Position = (1 - s - t)*p0 + s*p1 + t*p2;
  #elif GEOMETRY_ORDER == 2
  vec3 p3 = gl_in[3].gl_Position.xyz;
  vec3 p4 = gl_in[4].gl_Position.xyz;
  vec3 p5 = gl_in[5].gl_Position.xyz;

  float phi0 =  s*-3.0-t*3.0+s*t*4.0+(s*s)*2.0+(t*t)*2.0+1.0;
  float phi1 =  -s+(s*s)*2.0;
  float phi2 =  -t+(t*t)*2.0;
  float phi3 =  s*t*4.0;
  float phi4 =  t*(s+t-1.0)*-4.0;
  float phi5 =  -s*(s*4.0+t*4.0-4.0);

  v_Position = phi0 * p0 + phi1 * p1 + phi2*p2 + phi3*p3 + phi4*p4 + phi5*p5;
  #else
  #error "unsupported geometry order"
  #endif

  v_ParameterTess = vec2(s,t);
}

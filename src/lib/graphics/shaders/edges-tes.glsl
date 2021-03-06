//#version 410

layout (isolines, equal_spacing) in;

out vec3 v_Position;
out vec2 v_ParameterTess;

flat in float t_clip[];
flat out float v_clip;

void main() {

  float s = gl_TessCoord.x;
  float t = 1 - s;

  vec3 p0 = gl_in[0].gl_Position.xyz;
  vec3 p1 = gl_in[1].gl_Position.xyz;

  #if GEOMETRY_ORDER == 1
  v_Position = t*p0 + s*p1;
  #elif GEOMETRY_ORDER == 2

  vec3 p2 = gl_in[2].gl_Position.xyz;
  float phi0, phi1, phi2;
  phi0 =  s*-3.0+(s*s)*2.0+1.0;
  phi1 =  -s+(s*s)*2.0;
  phi2 =  s*4.0-(s*s)*4.0;
  v_Position = phi0 * p0 + phi1 * p1 + phi2*p2;

  #else
  #error "unsupported geometry order"
  #endif

  v_clip = t_clip[0];
  v_ParameterTess = vec2(s,t);
}

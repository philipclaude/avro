#version 410

layout (triangles, equal_spacing, ccw) in;

uniform mat4 u_ModelViewProjectionMatrix;

out vec3 v_Position;

void main() {

  float u = gl_TessCoord.x;
  float v = gl_TessCoord.y;

  vec3 P100 = gl_in[0].gl_Position.xyz;
  vec3 P010 = gl_in[1].gl_Position.xyz;
  vec3 P001 = gl_in[2].gl_Position.xyz;

  vec3 pos = (1 - u - v) * P100 + u*P010 + v*P001;


  gl_Position = u_ModelViewProjectionMatrix vec4(pos,1.0);
  v_Position  = pos;
}

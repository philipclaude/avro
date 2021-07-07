#version 410

layout (vertices = 3) out;

uniform int u_level;

void main() {
  gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;

  gl_TessLevelOuter[0] = u_level;
  gl_TessLevelOuter[1] = u_level;
  //gl_TessLevelInner[0] = u_level;
}

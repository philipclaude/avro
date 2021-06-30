#version 410

layout (vertices = 6) out;

void main() {
  gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;

  int level = 10;

  gl_TessLevelOuter[0] = level;
  gl_TessLevelOuter[1] = level;
  gl_TessLevelOuter[2] = level;

  gl_TessLevelInner[0] = level;
}

#version 410

#if GEOMETRY_ORDER == 1
#define NB_NODES 3
#elif GEOMETRY_ORDER == 2
#define NB_NODES 6
#elif GEOMETRY_ORDER == 3
#define NB_NODES 10
#elif GEOMETRY_ORDER == 4
#define NB_NODES 15
#else
#error "unsupported geometry order"
#endif

layout (vertices = NB_NODES) out;

uniform int u_level;

flat in float  v_clip[];
flat out float t_clip[];

void main() {
  gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;

  gl_TessLevelOuter[0] = u_level;
  gl_TessLevelOuter[1] = u_level;
  gl_TessLevelOuter[2] = u_level;

  gl_TessLevelInner[0] = u_level;

  t_clip[gl_InvocationID] = v_clip[gl_InvocationID];

}

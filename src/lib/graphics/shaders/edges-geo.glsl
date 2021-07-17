#version 410 core

layout (lines) in;
uniform mat4 u_ModelViewProjectionMatrix;

in vec3 v_Position[]; // receive the vertex coordinates from the vertex shader or tessellation shader

out vec2 v_Parameter; // pass the parameter coordinates to the fragment shader
out vec3 x_Position;

in vec2 v_ParameterTess[];

layout (line_strip , max_vertices = 2) out;

flat out float g_clip;
flat in float v_clip[];

uniform int u_clip;
uniform vec3 u_clip_center;
uniform vec3 u_clip_normal;

float
visible( in vec3 x ) {
  float result = 1.0;
  if (dot( x - u_clip_center , u_clip_normal ) < 0.0)
    result = -1.0;
  return result;
}

void main() {

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[0],1.0);
  v_Parameter = v_ParameterTess[0];
  x_Position  = v_Position[0];
  g_clip      = (u_clip == 2) ? v_clip[0] : visible( v_Position[0] );
  EmitVertex();

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[1],1.0);
  v_Parameter = v_ParameterTess[1];
  x_Position  = v_Position[1];
  g_clip      = (u_clip == 2) ? v_clip[1] : visible( v_Position[1] );

  // i'm not sure if this need to be set after every vertex, or just once before the end of the primitive
  // it at least needs to be set before the last EmitVertex()
  gl_PrimitiveID = gl_PrimitiveIDIn;

  EmitVertex();

  EndPrimitive();

}

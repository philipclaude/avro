#version 410 core

layout (lines) in;
uniform mat4 u_ModelViewProjectionMatrix;

in vec3 v_Position[]; // receive the vertex coordinates from the vertex shader or tessellation shader

out vec2 v_Parameter; // pass the parameter coordinates to the fragment shader

in vec2 v_ParameterTess[];

layout (line_strip , max_vertices = 2) out;

void main() {

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[0],1.0);
  v_Parameter = v_ParameterTess[0];
  EmitVertex();

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[1],1.0);
  v_Parameter = v_ParameterTess[1];

  // i'm not sure if this need to be set after every vertex, or just once before the end of the primitive
  // it at least needs to be set before the last EmitVertex()
  gl_PrimitiveID = gl_PrimitiveIDIn;

  EmitVertex();

  EndPrimitive();

}
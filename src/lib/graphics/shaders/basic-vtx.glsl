#version 330
layout (location = 0 ) in vec3 a_Position;

uniform mat4 u_ModelViewProjectionMatrix;

out vec3 v_Position;

void main()
{
  v_Position  = a_Position; // pass the position to the geometry shader
  gl_Position = u_ModelViewProjectionMatrix*vec4(a_Position,1.0);
}

#version 330 core

layout (triangles) in;
uniform mat4 u_ModelViewProjectionMatrix;

in vec3 v_Position[]; // receive the vertex coordinates from the vertex shader

out vec2 v_Parameter; // pass the parameter coordinates to the fragment shader

#if 1
layout (triangle_strip , max_vertices = 3) out;

void main() {

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[0],1.0);
  v_Parameter = vec2(0.,0.);
  EmitVertex();

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[1],1.0);
  v_Parameter = vec2(1.,0.);
  EmitVertex();

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[2],1.0);
  v_Parameter = vec2(1.,1.);
  EmitVertex();

  EndPrimitive();

}

#else

layout (line_strip , max_vertices = 4) out;

void main() {

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[0],1.0);
  v_Parameter = vec2(0.,0.);
  EmitVertex();

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[1],1.0);
  v_Parameter = vec2(1.,0.);
  EmitVertex();

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[2],1.0);
  v_Parameter = vec2(1.,1.);
  EmitVertex();

  gl_Position = u_ModelViewProjectionMatrix*vec4(v_Position[0],1.0);
  v_Parameter = vec2(0.,0.);
  EmitVertex();

  EndPrimitive();
  
}

#endif